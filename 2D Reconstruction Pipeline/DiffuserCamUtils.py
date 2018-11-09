
import numpy as np
import scipy.io as io
import scipy.sparse as sps
import matplotlib.pyplot as plt
from PIL import Image
from numpy.fft import fft2, ifft2, fftn, ifftn, fftshift, ifftshift
from numpy.linalg import norm
import math

"""
Utility functions for implementing the DiffuserCam code.
"""

"""
Options class for storing variables. Possible Fields (* parameters are required):

################################################################################
############################## GENERAL PARAMETERS ##############################
################################################################################

*type: 2d or 3d.
*algorithm: gd, nesterov, fista, or admm.
*iterations: max number of iterations.
gamma: step size per iteration.
eps: when cost < eps, solver terminates.
zsize: Size of zstack. If 2d, equals 1.
del_pixels: If true, will delete extra pixels from image. Default False.
psf: point spread function of system.
b: image taken.

################################################################################
##############################  CROP  PARAMETERS  ##############################
################################################################################

crop2d: 2d cropping function (crops all frames of a 3d stack).
crop3d: 3d cropping function used in forward model (performs crop2d then crops
        to the z = 0 frame).
pad2d: 2d padding function (pads all frames of a 3d stack).
pad3d: 3d padding function used in forward model (embeds pad2d into z = 0 frame
        and fills all other frames with zeros).
up_shape: unpadded shape of the 3d psf. If 2d problem, only one z dimension.
pad_shape: padded shape of the 3d psf.

################################################################################
############################  AUTOTUNE  PARAMETERS  ############################
################################################################################

autotune: true if autotuning during admm for mu paramters. Default True.
beta: threshold value for autotuning.
alpha: increment/decrement value for autotuning.

################################################################################
##############################  ADMM  PARAMETERS  ##############################
################################################################################
mu_i: mu parameters for admm.
tau: regularizer for TV norm.
"""
class Options:

    METHODS = ['gd', 'nesterov', 'fista', 'admm']

    def __init__(self, type, algorithm, iterations):
        assert type == '2d' or type == '3d'
        self.type = type
        assert algorithm in self.METHODS
        self.method = algorithm
        self.max_itr = iterations
        if type == '2d':
            self.zsize = 1
        self.eps = 0
        self.del_pixels = False
        self.autotune = True
        self.alpha = 0
        self.beta = 0

"""
Initialize the DiffuserCam.
"""
#def initialize(image_file, psf_file, f_lat = 8, f_ax = 4, type = 'pco', color = 'rgb', dim = '2d'):
def initialize(image_file, psf_file, f_lat = 1, f_ax = 1, type = 'pco', color = 'rgb', dim = '2d'):
    try:
        im_type = image_file[-3:]
        psf_type = psf_file[-3:]
        image = np.array(Image.open(image_file))
        psf = io.loadmat(psf_file)['psf'] if psf_type == 'mat' else np.array(Image.open(psf_file))
    except IOError as ex:
        print("I/O error: " + str(ex.strerror))
    if dim == '2d':         #embed into a 3d array
        psf = np.expand_dims(psf, 2)
    psf_bg = np.mean(image[0 : 15, 0 : 15]) #102
    image_bg = np.mean(image[0 : 15, 0 : 15])         #should be around 100

    psf_down = downsample_ax(psf - psf_bg, f_lat)
    image = downsample_ax(image - image_bg, f_lat)

    if dim == '3d':
        psf_down = downsample_lat(psf_down, f_ax)

    image /= np.max(image)
    psf_down /= norm(psf_down)

    return psf_down, image

"""
Get cropping and padding functions to use based on PSF shape. Padding will embed
the image into one with doubled dimensions. If N is provided, will zero out N
pixels.
"""
def get_crop_pad(psf, N = 0):

    up3d = np.array(psf.shape)
    up_shape = up3d[0 : 2]      #unpadded shape. don't care about zstack
    padded_shape = 2 * up_shape
    pad_shape_3d = (padded_shape[0], padded_shape[1], up3d[2])      #pad x and y, not z.

    start_crop = (padded_shape - up_shape) // 2
    end_crop = start_crop + up_shape

    # generate random incides to crop out
    (X, Y) = up_shape
    max_idx = X * Y
    rng = np.arange(max_idx)
    np.random.shuffle(rng)
    pixels = rng[:N]
    pixels = [(p % X, p // X) for p in pixels]

    # del_pixels: if true, will randomly delete N pixels.
    def crop2d(x, del_pixels = False):
        if len(x.shape) == 2:
            cropped = x[start_crop[0] : end_crop[0], start_crop[1] : end_crop[1]]
            if del_pixels and N > 0:
                cropped = pix_crop(cropped)
        elif len(x.shape) == 3:
            zsize = x.shape[2]
            cropped = np.zeros((up_shape[0], up_shape[1], zsize), 'float32')
            for z in range(zsize):
                cropped[:, :, z] = np.real(x[start_crop[0] : end_crop[0], start_crop[1] : end_crop[1], z])
                if del_pixels and N > 0:
                    cropped[:, :, z] = pix_crop(cropped[:, :, z])
        else:
            raise Exception('Object to crop must be 2d or 3d.')
        return cropped

    def pad2d(x, del_pixels = False):
        if del_pixels:
            x = pix_crop(x)
        if len(x.shape) == 2:
            padded = np.zeros(padded_shape, np.complex64)
            padded[start_crop[0] : end_crop[0], start_crop[1] : end_crop[1]] = x
        elif len(x.shape) == 3:
            padded = np.zeros((padded_shape[0], padded_shape[1], x.shape[2]), np.complex64)
            for z in range(x.shape[2]):
                padded[start_crop[0] : end_crop[0], start_crop[1] : end_crop[1], z] = x[:, :, z]
        else:
            raise Exception('Object to pad must be 2d or 3d.')
        return padded

    # 3d crop used in forward model. Performs a 2d crop over all slices and then
    # returns the z = 0 slice.
    def crop3d(x, del_pixels = False):
        return crop2d(x, del_pixels)[:, :, 0]

    # embeds 2d padded slice of x into the z = 0 slice of a 3d array with shape
    # pad_shape_3d
    def pad3d(x, del_pixels = False):
        padded = np.zeros(pad_shape_3d)
        padded[:, :, 0] = np.real(pad2d(x, del_pixels))         #embed the padded figure into the stack
        return padded

    def pix_crop(x):
        for (i, j) in pixels:
            x[i, j] = 0
        return x

    return crop2d, crop3d, pad2d, pad3d, pix_crop

"""
Returns all relevant operators for use in gradient descent.
h: padded psf.
crop: cropping operator.
pad: padding operator.
"""
def get_ops(h, crop2d, pad2d, crop3d, pad3d, up_shape):
    #for things which are not admm, put back in norm = 'ortho'
    H = fftn(ifftshift(h), norm = 'ortho')
    Hstar = np.conj(H)
    A = lambda x : np.real(crop3d(fftshift(ifftn(H * fftn(ifftshift(x), \
                        norm = 'ortho'), norm = 'ortho'))))
    AH = lambda x : np.real(fftshift(ifftn(Hstar * fftn(ifftshift(pad3d(x)), \
                        norm = 'ortho'), norm = 'ortho')))
    return A, AH

"""
Sets all negative values of an np array to 0. Used as a proximal for gradient descent.
"""
def non_negative(M):
    M[M < 0] = 0
    return M

"""
Minimize the cost function.

grad: operator to evaluate gradient at each iteration.
error: operator to evaluate error in objective function at each iteration.
prox: proximal to use at each iteration.
opt: options object with switches.
init: default initilization of object.
live_display: if True, will live update image of object every 10 iterations.
"""
def solver(grad, error, prox, opt, init = None, live_display = True):
    if init:
        if init.shape == opt.pad_shape:
            x0 = init
        elif init.shape == opt.up_shape:
            x0 = opt.pad(init)
        else:
            raise Exception('Shape of init matrix must be either padded or unpadded shape.')
    else:
        x0 = np.real(opt.pad2d(np.ones(opt.up_shape, dtype = 'float32')))   #initialize as a block of ones.

    #initialize live display
    if live_display:
        fig = plt.figure()
        ax = plt.gca()
        def update_display(i, x, mu_i = None):
            NUM_ITR = 20
            if i % NUM_ITR == 0:
                if i != 0:
                    disp_string = 'Iteration number ' + str(i)
                    if mu_i != None and opt.autotune:
                        for i in range(3):
                            disp_string += '. mu' + str(i + 1) + ': ' + \
                                                    '{:.1e}'.format(mu_i[i])
                    print(disp_string)
                x_display = np.mean(opt.crop2d(x), axis = 2)    #show avg. across z slices
                im = ax.imshow(x_display, cmap = 'gray')
                fig.canvas.draw()
        update_display(0, x0)
    else:
        update_display = lambda i, x: None

    if opt.method == 'gd':
        algorithm = gd
    elif opt.method == 'nesterov':
        algorithm = nesterov
    elif opt.method == 'fista':
        algorithm = fista
    elif opt.method == 'admm':
        algorithm = admm
    else:
        raise Exception('Invalid method parameter entered.')

    x, error_list = algorithm(grad, error, prox, opt, x0, update_display)
    return x, error_list

"""
Implementation of simple gradient descent.
"""
def gd(grad, error, prox, opt, x, update_display):
    i, error_list, e = 0, [], error(x)
    error_list.append(e)
    while i < opt.max_itr and (not opt.eps or e > opt.eps):
        y = x - opt.gamma * grad(x)
        x = prox(y)

        e = error(x)
        error_list.append(e)
        i += 1
        update_display(i, x)
    return x, error_list

"""
Implementation of Nesterov-accelerated gradient descent.
"""
def nesterov(grad, error, prox, opt, x, update_display):
    i, error_list, e = 0, [], error(x)
    error_list.append(e)
    p, mu = 0, .9
    while i < opt.max_itr and (not opt.eps or e > opt.eps):
        pk = p
        p = mu * p - opt.gamma * grad(x)
        y = x - mu * pk + (1 + mu) * p
        x = prox(y)

        e = error(x)
        error_list.append(e)
        i += 1
        update_display(i, x)
    return x, error_list

"""
Implementation of FISTA algorithm.
"""
def fista(grad, error, prox, opt, x, update_display):
    i, error_list, e = 0, [], error(x)
    error_list.append(e)
    xk_minus1 = x
    tk = 1
    while i < opt.max_itr and (not opt.eps or e > opt.eps):
        x -= opt.gamma * grad(x)
        xk = prox(x)
        tk1 = (1 + np.sqrt(1 + 4 * tk * tk)) / 2
        x = xk + (tk - 1) * (xk - xk_minus1) / tk1
        xk_minus1 = xk
        tk = tk1

        e = error(x)
        error_list.append(e)
        i += 1
        update_display(i, x)
    return x, error_list

"""
Implementation of ADMM. This implementation takes the crop from the aperture
into account and the non-negativity constraint, so we have added extra regularizers.
"""
def admm(grad, error, prox, opt, x, update_display):
    # initialize parameters
    h, b = opt.psf, opt.b
    zsize = opt.zsize           #size of zstack of psfs
    crop2d, pad2d, crop3d, pad3d = opt.crop2d, opt.pad2d, opt.crop3d, opt.pad3d
    gamma = opt.gamma   #step size
    size = opt.pad_shape
    del_pixels = opt.del_pixels     #if del_pixels, then C = crop3d + image crop
    autotune = opt.autotune
    alpha, beta = opt.alpha, opt.beta   #autotune parameters

    # get regularization parameters.
    if hasattr(opt, 'mu1'):
        mu1, mu2, mu3 = opt.mu1, opt.mu2, opt.mu3
        tau = opt.tau
    else:
        mu1, mu2, mu3 = 1e-4, 1e-3, 1e-2       # tune these. Working with 1e-4, 1e-4, tau = 1e-3
        tau = 1e-3                  #regularization parameter on TV norm
    print('mu1 = ' + str(mu1) + ', mu2 = ' + str(mu2) + ', mu3 = ' + str(mu3) + \
                    ', tau = ' + str(tau))
    if autotune:
        print('Autotuning with alpha = ' + str(alpha) + ', beta = ' + str(beta))

    # Flip h for convolution and construct convolution matrices.
    h = np.roll(np.flip(h, axis = 2), 1, axis = 2)
    H = fftn(ifftshift(h))
    HT = np.conj(H)
    DTD = computeDTD(size)
    HTH = HT * H
    J = mu1 * HTH + mu2 * DTD + mu3           #denom for x update

    #initialize matrices for crop.
    CTb = pad3d(b, del_pixels)
    CTC = pad3d(np.ones(b.shape, 'float32'), del_pixels)
    K = CTC + mu1           #denom for nu update

    #initialize M and MT. forward model is C * M(x)
    M = lambda x : np.real(fftshift(ifftn(H * fftn(ifftshift(x)))))
    MT = lambda x : np.real(fftshift(ifftn(HT * fftn(ifftshift(x)))))

    #initialize variables.
    x, xp = np.zeros(size, 'float32'), np.zeros(size, 'float32')
    xi, rho = np.zeros(size, 'float32'), np.zeros(size, 'float32')
    eta1, eta2, eta3 = np.zeros(size, 'float32'), np.zeros(size, 'float32'), np.zeros(size, 'float32')
    dxp, dyp, dzp = D(xp)
    Mxp = M(xp)

    #begin iteration
    i, error_list = 0, []
    e = error(xp)
    error_list.append(e)
    while i < opt.max_itr and (not opt.eps or e > opt.eps):

        #Store previous Mx for residual calculations
        Mx = Mxp
        dx, dy, dz = dxp, dyp, dzp

        # update u. Should be fine for 2d, because dz will be 0 --> z3 = 0
        z1, z2, z3 = dx + eta1 / mu2, dy + eta2 / mu2, dz + eta3 / mu2
        zmod = np.sqrt(z1 * z1 + z2 * z2 + z3 * z3)
        zmod[zmod <= 0] = 1       #don't divide by 0
        zmod = np.real(zmod)
        zmod = np.maximum(zmod - tau / mu2, np.zeros(size, 'float32')) / zmod
        u1, u2, u3 = z1 * zmod, z2 * zmod, z3 * zmod

        # update nu
        y = xi + mu1 * Mx + CTb
        nu = y / K

        #update w
        w = np.maximum(rho / mu3 + x, np.zeros(size, 'float32'))

        # update x
        r = DT(mu2 * u1 - eta1, mu2 * u2 - eta2, mu2 * u3 - eta3) + MT(mu1 * nu - xi) + (mu3 * w - rho)
        xp = np.real(fftshift(ifftn(fftn(ifftshift(r)) / J)))    #invert in Fourier space
        Mxp = M(xp)

        # update xi and mu1
        dxi = Mxp - nu
        xi = xi + mu1 * gamma * dxi
        if autotune:
            r1 = norm(dxi)
            s1 = mu1 * norm(Mxp - Mx)
            mu1, mu1_update = param_update(mu1, beta, alpha, r1, s1)

        #update eta and mu2
        dxp, dyp, dzp = D(xp)
        deta1, deta2, deta3 = dxp - u1, dyp - u2, dzp - u3
        eta1, eta2, eta3 = eta1 + mu2 * gamma * deta1, eta2 + mu2 * gamma * deta2, eta3 + mu2 * gamma * deta3
        if autotune:
            r2 = np.sqrt(norm(deta1) ** 2 + norm(deta2) ** 2 + norm(deta3) ** 2)
            s2 = mu2 * np.sqrt(norm(dxp - dx) ** 2 + norm(dyp - dy) ** 2 + norm(dzp - dz) ** 2)
            mu2, mu2_update = param_update(mu2, beta, alpha, r2, s2)

        #update rho and mu3
        dw = xp - w
        rho = rho + mu3 * gamma * dw
        if autotune:
            r3 = norm(dw)
            s3 = mu3 * norm(xp - x)
            mu3, mu3_update = param_update(mu3, beta, alpha, r3, s3)

            #if mus have been updated, update matrices as well.
            mu_update = mu1_update or mu2_update or mu3_update
            if mu_update:
                J = mu1 * HTH + mu2 * DTD + mu3
                K = CTC + mu1

        x = xp

        # append error and update display
        e = error(x)
        error_list.append(e)
        update_display(i, x, (mu1, mu2, mu3))
        i += 1
    return x, error_list

"""
Returns horizontal and vertical differences at each pixel. This is the finite
difference operator D_i x for use in ADMM.
x: unvectorized image to get differences from.
"""
def D(x):
    xroll = np.roll(x, -1, axis = 0)
    yroll = np.roll(x, -1, axis = 1)
    zroll = np.roll(x, -1, axis = 2)
    return xroll - x, yroll - x, zroll - x

"""
Adjoint of the difference operator hv_diff for use in ADMM.
dxi: Difference matrix along axis xi, for xi = x, y, or z
"""
def DT(dx, dy, dz):
    return (np.roll(dx, 1, axis = 0) - dx) + (np.roll(dy, 1, axis = 1) - dy) \
                    + (np.roll(dz, 1, axis = 2) - dz)

"""
Computes DT * D as a matrix multiplication, where D is the difference operator.
shape: of resulting matrix. Should be 3 dimensional.
"""
def computeDTD(shape):
    x = np.zeros(shape)
    if shape[2] == 1:           #then it is 2d
        x[0, 0, 0] = 4
        indices = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
    else:
        x[0, 0, 0] = 6
        indices = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    for (i, j, k) in indices:
        x[i, j, k] = -1
    DTD = np.abs(fftn(x))
    return DTD

"""
Update mu parameters in admm. Used in autotuning.
"""
def param_update(mu, beta, alpha, r, s):
    updated = False
    if r > beta * s:
        mu *= alpha
        updated = True
    elif r * beta < s:
        mu /= alpha
        updated = True
    return mu, updated

"""
Show PSF and input image.
"""
def show_data(psf, image):
    plt.figure()

    plt.subplot(121)
    plt.imshow(psf, cmap = 'gray')
    plt.title('PSF')

    plt.subplot(122)
    plt.imshow(image, cmap = 'gray')
    plt.title('Image Data')

    plt.show()

"""
Computes the gradient of cost function f(x) = 1/2 ||Ax - b||^2.
x: input image.
A: transfer function for forward model.
AH: adjoint of A.
b: measured image.
"""
def grad(x, A, AH, b):
    return AH(A(x) - b)

"""
Computes the cost function ||b - A(x)||_2 + tau * |Dx|_1 at each step.
"""
def objective(x, A, b, tau):
    return norm(A(x) - b) + tau * np.sum(np.abs(D(x)))

"""
Downsample image by given factor along the transverse direction.
"""
def downsample_ax(img, factor):
    n = int(np.log2(factor))
    for i in range(n):
        if len(img.shape) == 2:
            img = .25 * (img[::2, ::2] + img[1::2, ::2]
                + img[::2, 1::2] + img[1::2, 1::2])
        else:
            img = .25 * (img[::2, ::2, :] + img[1::2, ::2, :]
                + img[::2, 1::2, :] + img[1::2, 1::2, :])
    return img

"""
Downsample a zstacked image by a given factor in the z direction.
"""
def downsample_lat(img, factor):
    n = int(np.log2(factor))
    for i in range(n):
        img = .5 * (img[:, :, ::2] + img[:, :, 1::2])
    return img
