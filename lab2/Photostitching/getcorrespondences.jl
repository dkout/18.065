using GeometryTypes, Colors, GLWindow, GLAbstraction
using FileIO, Reactive, GLVisualize, GLFW
import GLVisualize: mm, layoutscreens, button, IRect
import GLAbstraction: singlepressed
FocusWindow(window::GLFW.Window) = ccall( (:glfwFocusWindow, GLFW.lib), Void, (GLFW.WindowHandle,), window)

function getcorrespondences(img1, img2, n)
    mw, mh = GLWindow.primarymonitorresolution()
    mw, mh = max(mw - 50, 200), max(mh - 50, 200) # don't fill complete screen
    imw, imh = size(img1, 2) + size(img2, 2), max(size(img1, 1), size(img2, 1))
    wres, hres = min(mw, imw), min(mh, imh)
    window = glscreen("Pick Points", resolution = (wres, hres))
    @async GLWindow.waiting_renderloop(window)
    nimages = 2
    layout = [
        [Signal("Image $(nimages - i + 1)") for i = 1:nimages],
    ]

    screens = layoutscreens(
        window, layout,
        stroke = (1mm, RGBA(0.3f0, 0.3f0, 0.3f0)),
        text_color = RGBA(0.2f0,0.2f0,0.2f0),
        title_height = 8mm
    )

    image_screens = screens[1]
    img_pictures = [img1, img2]
    startpoints = [Signal(Point2f0[]) for i = 1:nimages]

    img_transform = Signal{SimpleRectangle{Int}}[]
    for (screen, im) in zip(image_screens, img_pictures)
        prim = map(screen.area) do a
            hi, wi = size(im)
            w1, w2 = Vec(wi, hi), Vec(a.w, a.h)
            ratio = w2 ./ w1
            if ratio[1] < ratio[2]
                IRect(0, 0, a.w, hi * ratio[1])
            else
                IRect(0, 0, wi * ratio[2], a.h)
            end
        end
        push!(img_transform, prim)
        _view(visualize(
            im, primitive = prim
        ), screen, camera = :fixed_pixel)
    end
    r1, r2 = map(value, img_transform)
    resize!(window, r1.w + r1.w, max(r1.h, r2.h))
    picking_screen = Signal(1)
    preserve(map(picking_screen) do i
        if i != 0
            for j = 1:nimages
                message = ""
                if i == j
                    image_screens[j].stroke = (2mm, RGBA(1.0f0, 0.3f0, 0.3f0))
                    selected = length(value(startpoints[i]))
                    if i == 1
                        message = "Select a point ($(n-selected) left)"
                    else
                        message = "Select a corresponding point ($(n-selected) left)"
                    end
                else
                    image_screens[j].stroke = (1mm, RGBA(0.3f0, 0.3f0, 0.3f0))
                    message = "Image $(nimages - i + 1)"
                end
                push!(layout[1][j], message)
            end
        end
        nothing
    end)
    for (i, screen) in enumerate(image_screens)
        showcross = map(screen.inputs[:mouseinside], picking_screen) do minside, pidx
            minside && pidx == i
        end
        m = 5000; g = 1mm
        buff = fill(Point2f0(0), 8)
        positions = foldp(buff, screen.inputs[:mouseposition]) do v0, mp
                v0[1] = (-m, mp[2]); v0[2] = (mp[1] - g, mp[2])
                v0[3] = (mp[1] + g, mp[2]); v0[4] = (m, mp[2])
                v0[5] = (mp[1], -m); v0[6] = (mp[1], mp[2] - g)
                v0[7] = (mp[1], mp[2] + g); v0[8] = (mp[1], m)
                v0
        end
        _view(visualize(
            positions, :linesegment, visible = showcross
        ), screen, camera = :fixed_pixel)
    end
    # Visualize the points on each screen
    for (i, points) in enumerate(startpoints)
        _view(visualize(
            ('+', points), # points centered, 10pixel radius
            scale = Vec2f0(5mm), offset = Vec2f0(-2.5mm),
            color = RGBA(1f0, 0.3f0, 0.3f0),
        ), image_screens[i], camera = :fixed_pixel)
    end

    @materialize mouse_buttons_pressed, mouseinside = window.inputs

    preserve(map(mouse_buttons_pressed) do mbp
        idx = value(picking_screen)
        if idx != 0
            screen = image_screens[idx]
            if singlepressed(mbp, GLFW.MOUSE_BUTTON_LEFT) && value(screen.inputs[:mouseinside])
                mpos = value(screen.inputs[:mouseposition])
                points_s = startpoints[idx]
                push!(points_s, push!(value(points_s), mpos))
                push!(picking_screen, mod1(idx + 1, nimages))
            end
        end
        return
    end)

    # transform all points back to their pixel space
    points_p = map(zip(startpoints, img_transform, img_pictures)) do pointstransformareaim
        ps, rect, im = pointstransformareaim
        map(ps, rect) do ps, a
            hi, wi = size(im)
            w1, w2 = Vec(wi, hi), Vec(a.w, a.h)
            scale = w1 ./ w2
            map(ps) do p
                p = scale .* Point2f0(p[2], p[1])
                (p[2], hi - p[1])
            end
        end
    end
    println("start picking points: ")
    img = zeros(RGB{N0f8}, 10, 10)
    FocusWindow(GLWindow.nativewindow(window))
    while isopen(window)
        tic()
        if length(value(last(points_p))) == n
            img = RGB{N0f8}.(GLWindow.screenbuffer(window))
            println("done")
            destroy!(window)
            break
        end

        yield()
        t = toq()
        # try to go for 60 fps
        # sleep the remaining time
        t = max((1/60) - t, 0)
        t >= 0.001 && sleep(t) # 0.001 == minimum time sleep can do
    end
    points = map(points_p) do points
        ps = value(points)
        [ps[i][j] for i = 1:length(ps), j = 1:2]
    end
    (points[1], points[2]), img
end


# names = ["law"] ## sample file -- tag the logos in the window in first image
# im_name = names[1]
# # test with:
# inpath1 = "$(im_name)1.jpg"
# inpath2 = "$(im_name)2.jpg"

# ## Load images and convert to an Array
# im1 = convert(Array, load(inpath1))
# im2 = convert(Array, load(inpath2))
# getcorrespondences(im1, im2, 5)
