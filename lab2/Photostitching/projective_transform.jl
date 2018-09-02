function projective_transform(XY, XYtil)
    #
    # Syntax: H = projective transform(XY, XYtil)
    #
    # Inputs: XY and XYtil are n x 2 matrices containing (x, y) coordinates
    # for n corresponding points in two coordinate systems
    #
    # Outputs: H is the unique 3 x 3 projective transformation matrix that maps
    # XY to XYtil with H[3, 3] = 1. That is, in the ideal case, the
    # following relationship should hold:
    #
    # tmp = [XY ones(n)] * H
    # XYtil = tmp[:, 1:2] ./ tmp[:, 3]

    n = size(XY, 1)
    # Construct A matrix
    A = zeros(2n, 9)
    for i in 1:n
        xy1i = [XY[i, :]; 1]
        A[2i - 1, :] = [zeros(3); -xy1i; XYtil[i, 2] * xy1i]
        A[2i, :] = [xy1i; zeros(3); -XYtil[i, 1] * xy1i]
    end

    v = svd(A, thin=false)[3][:, 9]
    H = reshape(v / v[9], 3, 3)
    return H
end
