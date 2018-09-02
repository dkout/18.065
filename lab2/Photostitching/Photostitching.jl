
# module names should be written in camel case (every first letter of a word should be uppercase)
# this will isolate this code
module Photostitching

using Interpolations

type ImageStitcher
    image
    tform
    xlim
    ylim
end

ImageStitcher(image, tform) = ImageStitcher(image, tform, Float64[], Float64[])

# in a module we need to export the names, to be imported via the `using` command
export ImageStitcher


"""
 Syntax:       Is, alpha = stitchImages(It);
               Is, alpha = stitchImages(...,'dim',dim,...);
               Is, alpha = stitchImages(...,'b0',b0,...);
               Is, alpha = stitchImages(...,'view',view,...);
               Is, alpha = stitchImages(...,'order',order,...);
               Is, alpha = stitchImages(...,'mode',mode,...);

 Inputs:       It is an n x 1 array of structs, where

                   It(i).image is an Nyi x Nxi x Nc uint8 image

                   It(i).tform is the 3 x 3 transformation matrix that
                   maps It(i).image into the coordinates specified by the
                   mode argument below

               [OPTIONAL] dim = [h, w] are the desired height and width,
               respectively, in pixels, of the stitched image. By default
               dim is set to preserve the input resolutions

               [OPTIONAL] b0 is a Nc x 1 vector of double values in
               [0, 255] specifying the background (i.e., missing data)
               color. The default value is b0 = 255 * ones(Nc,1)

               [OPTIONAL] view = {'left','center','right','default'}
               specifies what view to construct the stitched image from.
               The default is view = 'default' (world coordinates)

               [OPTIONAL] order = {'natural','reverse'} specifies the
               order in which to overlay the stitched images. In natural
               ordering, the closest image is displayed on top. The
               default is order = 'natural'

               [OPTIONAL] mode = {'absolute','relative'} specifies the
               type of transform used. When mode = 'absolute', It(i).tform
               should map to "world" coordinates. When mode = 'relative',
               It(i).tform should map It(i).image into the coordinates of
               It(i-1).image. Here, It(1).tform can be eye(3) or another
               transformation that maps to the desired world coordinates.
               The default value is mode = 'relative'

 Outputs:      Is is a h x w x Nc uint8 matrix containing the stitched
               image

               alpha is a h x w uint8 matrix containing the alpha channel
               (transparency) values for generating a nice transparent png
               of the stitched image

 Author:       Brian Moore
               brimoor@umich.edu

 Date:         September 30, 2015
"""
function stitchImages(It;
        view = "default",
        order = "natural",
        mode = "relative",
        dim = [],
        b0 = one(eltype(It[1].image))
    )

    # Convert to absolute coordinates, if necessary
    if mode == "relative"
        toAbsoluteCoordinates(It)
    end

    # Apply view
    applyView(It, view)

    # Apply order
    applyOrder(It, order)

    # Compute stitched image limits
    xlim, ylim = computeStitchedLimits(It)

    # Get sample points
    x, y, w, h = getSamplePoints(xlim, ylim, dim)

    # Stitch images
    Is = fill(RGBA(NaN, NaN, NaN, NaN), h, w)
    for i = length(It):-1:1
        overlayImage(Is, It[i], x, y)
    end

    # Fill background
    fillBackground(Is, b0)

    Is
end

export stitchImages

#--------------------------------------------------------------------------
function toAbsoluteCoordinates(It)
    # Map all transforms to coordinates of It(1)
    for i = 2:length(It)
        It[i].tform = It[i].tform * It[i - 1].tform;
    end
end

#--------------------------------------------------------------------------
function applyView(It, view)
    # Get transformed limits
    Ns = length(It)
    xc = zeros(Ns)
    for i = 1:Ns
        xlimi, _ = getOutputLimits(It[i].image, It[i].tform)
        xc[i] = mean(xlimi)
    end

    # Get ordering
    view = lowercase(view)
    if "left" == view
        # Left view
        idx = sort(xc)
    elseif "center" == view
        # Center view
        idx = sort(abs(xc - mean(xc)));
    elseif "right" == view
        # Right view
        idx = sort(xc, rev = true)
    elseif "default" == view
        # Use input ordering
        idx = 1:Ns
    else
        # Unsupported view
        error("Unsupported view $view")
    end

    # Apply ordering
    It[1:end] = It[idx]
    if view != "default"
        H1 = It[1].tform
        for i = 1:Ns
            It[i].tform = It[i].tform / H1
        end
    end
end
#--------------------------------------------------------------------------
function applyOrder(It, order)
    # Apply order
    order = lowercase(order)
    if order == "natural"
        # Natural ordering
        # Empty
    elseif order == "reverse"
        # Reverse ordering
        reverse!(It)
    else
        # Unsupported order
        error("Unsupported order: $order")
    end
end
#--------------------------------------------------------------------------
function computeStitchedLimits(It)
    # Compute limits
    minx = Inf;
    maxx = -Inf;
    miny = Inf;
    maxy = -Inf;
    for i = 1:length(It)
        xlimi, ylimi = getOutputLimits(It[i].image, It[i].tform)
        It[i].xlim = xlimi
        It[i].ylim = ylimi
        minx = min(minx, xlimi[1])
        maxx = max(maxx, xlimi[2])
        miny = min(miny, ylimi[1])
        maxy = max(maxy, ylimi[2])
    end
    xlim = [floor(Int, minx), ceil(Int, maxx)]
    ylim = [floor(Int, miny), ceil(Int, maxy)]
    xlim, ylim
end
#--------------------------------------------------------------------------
function getOutputLimits(I, H)
    # Compute limits of transformed image
    Ny, Nx = size(I)
    X = [1 Nx Nx 1]
    Y = [1 1 Ny Ny]
    Xt, Yt = applyTransform(X, Y, H)
    xlim = [minimum(Xt), maximum(Xt)]
    ylim = [minimum(Yt), maximum(Yt)]
    xlim, ylim
end


#--------------------------------------------------------------------------
function applyTransform(X, Y, H)
    # Apply transformation
    sz = size(X)
    n = length(X)
    tmp = [X[:] Y[:] ones(n)] * H
    Xt = reshape(tmp[:, 1] ./ tmp[:, 3], sz)
    Yt = reshape(tmp[:, 2] ./ tmp[:, 3], sz)
    Xt, Yt
end

#--------------------------------------------------------------------------
function applyInverseTransform(Xt, Yt, H)
    # Apply inverse transformation
    sz = size(Xt)
    n = length(Xt)
    tmp = [Xt[:] Yt[:] ones(n)] / H
    X = reshape(tmp[:, 1] ./ tmp[:, 3], sz)
    Y = reshape(tmp[:, 2] ./ tmp[:, 3], sz)
    X, Y
end
#--------------------------------------------------------------------------
function getSamplePoints(xlim, ylim, dim)
    # Get sample dimensions
    if isempty(dim)
        w = (diff(xlim) + 1)[1]
        h = (diff(ylim) + 1)[1]
    else
        h, w = dim
    end

    # Limit resolution to a reasonable value, if necessary
    MAX_PIXELS = 2000 * 2000
    w, h = limitRes(w, h, MAX_PIXELS)

    # Compute sample points
    x = linspace(xlim[1], xlim[2], w)
    y = linspace(ylim[1], ylim[2], h)
    x, y, w, h
end

#--------------------------------------------------------------------------
function limitRes(w, h, lim)
    if w * h <= lim
        # No rescaling needed
        return w, h
    end

    # Rescale to meet limit
    kappa = w / h
    w = round(Int, sqrt(lim * kappa))
    h = round(Int, sqrt(lim / kappa))
    warn("Output resolution too large, rescaling to $h x $w")
    w, h
end

#--------------------------------------------------------------------------
function overlayImage(Is, It, x, y)
    # Overlay image
    If = fillImage(It, x, y)
    mask = ~any(isnan.(If), 3)
    Isj = Is[:, :]
    Ifj = If[:, :]
    Isj[mask] = Ifj[mask]
    Is[:, :] = Isj
end

#--------------------------------------------------------------------------
function fillImage(It, x, y)
    # Parse inputs
    w = length(x)
    h = length(y)

    # Get active coordinates
    xIdx1 = findlast( x-> x <= It.xlim[1], x)
    xIdx2 = findfirst(x-> x >= It.xlim[2], x)

    yIdx1 = findlast( y-> y <= It.ylim[1], y)
    yIdx2 = findfirst(y-> y >= It.ylim[2], y)
    wa = xIdx2 + 1 - xIdx1
    ha = yIdx2 + 1 - yIdx1

    # Compute inverse transformed coordinates
    # part overlapping with It.image
    Xta = [x[i] for j in yIdx1:yIdx2, i in xIdx1:xIdx2]
    Yta = [y[j] for j in yIdx1:yIdx2, i in xIdx1:xIdx2]
    # coordinates of target image transformed to It.image
    Xa, Ya = applyInverseTransform(Xta, Yta, It.tform)

    # Compute active image
    Ia = fill(RGBA(NaN, NaN, NaN, NaN), ha, wa)
    imgip = interpolate(It.image, BSpline(Linear()), OnGrid())
    for xi = 1:wa, yi = 1:ha
        if xi <= size(Xa, 2) && yi <= size(Xa, 1)
            xt, yt = Xa[yi, xi], Ya[yi, xi]
            if xt <= size(imgip, 2) && yt <= size(imgip, 1) && yt >= 1 && xt >= 1
                Ia[yi, xi] = imgip[yt, xt]
            end
        end
    end
    
    # Embed into full image
    If = fill(RGBA(NaN, NaN, NaN, NaN), h, w)
    If[yIdx1:yIdx2, xIdx1:xIdx2] = Ia
    If
end
#--------------------------------------------------------------------------
function fillBackground(Is, b0)
    # Fill background
    mask = any(isnan.(Is), 3)
    Isj = Is[:, :]
    Isj[mask] = b0
    Is[:, :] = Isj

    # alpha = zeros(N0f8, size(mask))
    # alpha[~mask] = 1.0
end

# include GLVisualize point picking UI
include("getcorrespondences.jl")
export getcorrespondences

end

# using the above module and making the exported functions available
# The dot will search for the module in the current scope, as opposed to global scope
using .Photostitching
