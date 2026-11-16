# YARRAY: seismic-array processing

`shakelab.exploration.yarray` measures the waves crossing a seismic
array by frequency-wavenumber (f-k) analysis: phase velocity and
propagation azimuth as functions of frequency and time, and the signed
ellipticity of Rayleigh waves. It also tells, before or after a
survey, what an array can measure: its geometric response, the
statistical performance of the estimators and the frequency-velocity
domain where the analysis is reliable. The frequency domain
decomposition (FDD, also enhanced) uses the same cross-spectral
matrices to identify resonances, their shapes and their damping.

The same components are used in three ways:

- as a **library** from Python (`ArrayData`, `Beamformer`, `FDD`,
  results, statistics, `ArrayResponse` for planning);
- from the **command line** (`slab-yarray`, beamforming and array
  response);
- in **real time**, as a function of the SeedLink client
  (`BeamformingMonitor`, `FDDMonitor`).

This guide gives the theory (Part I), the use (Part II) and the
design, validation and limits (Part III). An overview of the package,
with its main uses and modules, is in the
[package README](../../shakelab/exploration/yarray/README.md); scripts
illustrating every way of use are in
[`examples/exploration/yarray/`](../../examples/README.md#exploration).
The guide derives from the four documents of the YARRAY development code
(August 2026: theoretical background, software architecture, user
manual, array response theory), updated to the implementation in
ShakeLab.

## Contents

Part I. Theory

1. [Fourier coefficients and cross-spectral matrices](#1-fourier-coefficients-and-cross-spectral-matrices)
2. [Plane waves and steering vectors](#2-plane-waves-and-steering-vectors)
3. [Estimators](#3-estimators)
4. [Components and Rayleigh three-component beamforming](#4-components-and-rayleigh-three-component-beamforming)
5. [Uncertainty of the power](#5-uncertainty-of-the-power)
6. [Spectral maxima and their density](#6-spectral-maxima-and-their-density)
7. [Geometric array response](#7-geometric-array-response)
8. [Statistical performance and operational limits](#8-statistical-performance-and-operational-limits)
9. [Frequency domain decomposition](#9-frequency-domain-decomposition)

Part II. Use

10. [Data](#10-data)
    - [Array file](#array-file)
11. [Configuration](#11-configuration)
12. [Processing and results](#12-processing-and-results)
    - [Frequency domain decomposition](#frequency-domain-decomposition)
    - [Result files](#result-files)
13. [Command line](#13-command-line)
14. [Real time](#14-real-time)
15. [Survey planning](#15-survey-planning)
16. [Practical guidance and processing times](#16-practical-guidance-and-processing-times)

Part III. Design, validation and limits

17. [Design choices](#17-design-choices)
18. [Validation](#18-validation)
19. [Current limitations and planned extensions](#19-current-limitations-and-planned-extensions)
20. [References](#20-references)

# Part I. Theory

The wavefield is assumed locally plane over the array and propagating
horizontally. Coordinates are east ($`x`$), north ($`y`$) and up, in
metres (or in any distance unit, consistently with velocities and
wavenumbers); azimuths are clockwise from north and give the direction
towards which the wave travels.

## 1. Fourier coefficients and cross-spectral matrices

Each record is cut into elementary **windows**, tapered (Tukey) and
transformed at the analysis frequencies. The windows have a fixed
duration, or a duration of $`n_c`$ cycles of each frequency,
$`T_w = n_c / f`$; they may overlap. Consecutive windows form
**covariance blocks**, which may also overlap: each block gives one
estimate of the cross-spectral matrix (CSM) and therefore one
spectrum, so that a long record gives a time series of results.

For window $`m`$ the coefficients of the $`N`$ channels at frequency
$`f`$ form a vector $`𝐱_m(f)`$. The CSM of a block of $`L`$
windows is

```math
𝐂(f) = \frac{1}{L} \sum_{m=1}^{L}
    𝐱_m(f)\, 𝐱_m^H(f)
```

The coefficient of each channel is scaled as a one-sided power
spectral density,

```math
x(f) = \sqrt{\frac{2 \Delta}{\sum_n w_n^2}}
    \sum_n w_n\, x_n\, e^{-2 \pi i f n \Delta}
```

($`w_n`$ taper, $`x_n`$ samples, $`\Delta`$ sampling interval), so
that the diagonal of $`𝐂`$ is the PSD of each channel and
Bartlett and Capon give the density of a plane wave in the units of
the data squared per Hz. Data not corrected for the instrument work
the same way: the power is then in counts squared per Hz, or
`(unit)**2/Hz` when the units are unknown.

Two transforms give the coefficients
(`shakelab.signals.crossspectral.coefficients`):

- **direct**: the transform evaluated at each requested frequency; any
  frequency, fixed or cycle windows;
- **RFFT**: one real FFT per window (fixed windows only), at the
  **nearest** bin (repeated bins are processed once, so the
  frequencies are the bins) or **interpolated** at the requested
  frequencies (zero padding and local cubic interpolation of the
  complex coefficients). Zero padding densifies the sampled transform
  but does not improve the resolution of the window.

**Effective number of windows.** Overlapping tapered windows are
correlated, so $`L`$ windows carry the information of fewer
independent ones (Percival and Walden, 1993, Section 6.17):

```math
L_{\mathrm{eff}} = \frac{L}
    {1 + 2 \sum_{q} \left(1 - \frac{q}{L}\right) \rho_q^2}
```

with $`\rho_q`$ the overlap of the taper with itself shifted by $`q`$
windows. $`L_{\mathrm{eff}}`$ is used by the standard errors of the
power.

## 2. Plane waves and steering vectors

A plane wave of angular wavenumber $`k`$ and azimuth $`\theta`$ has
the wavevector $`𝐤 = (k_x, k_y)`$,

```math
k_x = k \sin\theta, \qquad k_y = k \cos\theta, \qquad
k = \frac{2 \pi f}{c} = 2 \pi f s
```

($`c`$ phase velocity, $`s`$ slowness), and at station $`i`$, at
$`(x_i, y_i)`$, the phase factor

```math
a_i(𝐤) = e^{-i (k_x x_i + k_y y_i)}
```

The steering vectors $`𝐚`$ are evaluated on a polar **search
grid**, whose radial coordinate is velocity, slowness or wavenumber,
sampled linearly or logarithmically, with uniform or explicit
azimuths.

## 3. Estimators

**Bartlett** (conventional beamforming):

```math
P_B = \frac{𝐚^H 𝐂\, 𝐚}{N^2}
```

robust and inexpensive; its resolution is that of the array response
(Section 7).

**Capon** (minimum variance, Capon 1969), with trace-scaled diagonal
loading $`\lambda`$ (`diagonal_loading`) for numerical stability:

```math
𝐑 = 𝐂
    + \lambda\, \frac{\mathrm{tr}(𝐂)}{N}\, 𝐈,
\qquad
P_C = \frac{1}{𝐚^H 𝐑^{-1} 𝐚}
```

It suppresses the power not coherent with the tested wave and resolves
closer waves than Bartlett. A CSM of $`L`$ windows underestimates the
Capon power by about $`(L - N + 1) / L`$ (Capon and Goodman, 1970):
the bias does not move the maxima. With `strict=True` (default) Capon
and MUSIC require more windows per block than the dimension of the
matrix.

**MUSIC** (Schmidt, 1986): with $`𝐄_n`$ the eigenvectors of
the noise subspace (all but `music_sources`),

```math
P_M = \frac{1}{𝐚^H 𝐄_n 𝐄_n^H 𝐚}
```

a pseudo-spectrum without physical units: very sharp when the
subspace model holds, sensitive to the assumed number of sources.

## 4. Components and Rayleigh three-component beamforming

For each tested azimuth $`\theta`$ the horizontal motion (east $`E`$,
north $`N`$) is projected on the radial and transverse directions,

```math
R = E \sin\theta + N \cos\theta, \qquad
T = E \cos\theta - N \sin\theta
```

(transverse: radial rotated by +90 degrees, clockwise). The component
modes are: `vertical`, `radial`, `transverse`,
`three_component_separate` (the three spectra on the same grid),
`rayleigh_joint`, `rayleigh_angular` and `rayleigh_direct`.

**Rayleigh three-component beamforming** (RTBF, Wathelet et al.,
2018) combines, for each azimuth, the radial and vertical coefficients
in one matrix of dimension $`2N`$, steered by

```math
𝐬(e) =
    \begin{bmatrix} -i\, e\, 𝐚 \\ 𝐚 \end{bmatrix}
```

where $`e`$ is the signed ellipticity (radial over vertical amplitude,
positive for radial motion $`e \sin(\omega t)`$ with vertical
$`\cos(\omega t)`$). The radial component is positive in the direction
of propagation and the vertical one upwards, so that **$`e > 0`$ is
prograde and $`e < 0`$ retrograde** particle motion, the convention of
Wathelet et al. (2018) (other authors use the opposite sign). The
fundamental mode of a homogeneous half-space, retrograde at the
surface, gives $`e = -0.68`$ for a Poisson solid; for a soft layer
over a stiff half-space the fundamental mode is retrograde below the
peak of the H/V ratio, prograde between the peak and the following
trough, and retrograde again above it. Data whose vertical component
is positive downwards (or whose horizontal components are reversed)
give the opposite sign: check the polarity of the sensors. The Capon
power $`P`$ of $`𝐬(e)`$ is computed for every trial
ellipticity and the objective $`e^2 P^2`$ is maximized: the spectrum
is that maximum, and the ellipticity at the maximum is kept with each
spectral maximum. The ellipticity is sampled uniformly in its angle
$`\chi = \arctan e`$ strictly inside (-90, 90) degrees (default -85
to 85, 61 values), so that large ellipticities are reached without an
arbitrary bound on $`e`$; zero is excluded. Its power has the square
of the units of the power.

Wathelet (2024) shows that the 2018 formulation is distorted by
incoherent noise for vanishing or very large ellipticities and
proposes a normalized steering in the ellipticity angle
$`\xi = \arctan e`$, with the power itself maximized over $`\xi`$:

- **ARTBF** (`rayleigh_angular`): the radial-vertical matrix of
  dimension $`2N`$ of each azimuth, steered by $`𝐬_A(\xi)`$;
- **ARDS** (`rayleigh_direct`): one matrix of dimension $`3N`$ of the
  east, north and vertical motion per block (inverted once for
  Capon), steered by $`𝐬_D(\xi)`$:

```math
𝐬_A(\xi) =
    \begin{bmatrix} -i \sin\xi\, 𝐚 \\
                     \cos\xi\, 𝐚 \end{bmatrix},
\qquad
𝐬_D(\xi) =
    \begin{bmatrix} -i \sin\xi \sin\theta\, 𝐚 \\
                     -i \sin\xi \cos\theta\, 𝐚 \\
                     \cos\xi\, 𝐚 \end{bmatrix}
```

Both work with Capon and Bartlett. The steering vector is
$`𝐬 = 𝐀 𝐮`$ with
$`𝐮 = (\sin\xi, \cos\xi)^T`$ real and of unit length, so
that

```math
𝐬^H 𝐖 𝐬
    = 𝐮^T \mathrm{Re}(𝐀^H 𝐖 𝐀)
      \, 𝐮
    = m + d \cos(2 \xi - \varphi)
```

with $`𝐖`$ the inverse matrix (Capon) or the matrix
(Bartlett) and $`m`$, $`d`$, $`\varphi`$ given by the three entries of
the real 2 x 2 matrix
$`𝐐 = \mathrm{Re}(𝐀^H 𝐖 𝐀)`$.
The maximum over $`\xi`$ is therefore found exactly, by the
eigenvector of $`𝐐`$ (Capon: $`P = 1 / (m - d)`$; Bartlett:
$`P = (m + d) / N^2`$), with no grid of ellipticities and without the
Newton search proposed by Wathelet (2024). The angle covers
(-90, 90] degrees: purely vertical ($`\xi = 0`$) and purely
horizontal motion ($`\xi = 90`$ degrees, no finite $`e`$, stored as
None) are reached. The power is the density of the radial and
vertical motion together, $`|R|^2 + |Z|^2`$, in the units of the
power. The joint processing uses the same quadratic forms on its
ellipticity grid (identical results, more than ten times faster than
solving one system per ellipticity).

Capon underestimates the power by about $`(L - n + 1) / L`$ for a
matrix of dimension $`n`$ and $`L`$ windows (Section 19): more for the
matrices of dimension $`2N`$ and $`3N`$ than for single components,
so compare powers of different modes with Bartlett or with many
windows. The correction of the noise difference between vertical and
horizontal components (Wathelet, 2024, Section 3.6) is not
implemented yet (Section 19).

## 5. Uncertainty of the power

`statistical_uncertainty="empirical"` gives the standard error of the
power at every grid point (single components):

- Bartlett: the power of each window is computed separately and the
  variance of their mean is divided by $`L_{\mathrm{eff}}`$;
- Capon: first-order propagation (delta method) of the variance of the
  windows through
  $`\delta P = P^2\, 𝐛^H \delta𝐑\, 𝐛`$,
  $`𝐛 = 𝐑^{-1} 𝐚`$.

`geometric_uncertainty="expected_coherence"` (Bartlett) gives the
power expected over independent Gaussian errors of the station
coordinates, with covariance $`𝐒_i`$: the steering
coefficients are damped by the expected coherence
$`g_i = \exp(-𝐤^T 𝐒_i 𝐤 / 2)`$ and the
diagonal is kept exactly,

```math
P = \frac{(𝐚 \circ 𝐠)^H\, 𝐂\,
          (𝐚 \circ 𝐠)}{N^2}
    + \frac{1}{N^2} \sum_i \left(1 - g_i^2\right) C_{ii}
```

($`\circ`$ element-wise product). These are uncertainties of the
power, not of the position of a maximum.

## 6. Spectral maxima and their density

At each frequency, block and component the local maxima of the
spectrum (azimuth x radial neighbourhood, default 3 x 3, optional
minimum power relative to the largest) are ranked by power; rank 1 is
the strongest, not a mode. Each maximum keeps frequency, block and its
time (centre of the block), wavenumber and azimuth (velocity,
slowness, $`k_x`$, $`k_y`$ derived), power, its standard error and,
for the Rayleigh processing, the ellipticity.

The maxima of many blocks (a long record, or the analyses of the
monitor joined) scatter around the dispersion curves. Their density at
each frequency is estimated as a weighted mixture of Gaussian kernels
in wavenumber, log-slowness, log-velocity or, for the Rayleigh
processing, ellipticity angle $`\arctan e`$:

```math
D(x \mid f) = \sum_i w_i\, \mathcal{N}(x;\, x_i, \sigma_i)
```

with uniform, power or relative-power weights $`w_i`$ and kernels
as wide as the resolution of the search grid (default), Silverman's
rule or a fixed width. Its summary gives mode, median and a quantile
interval; its local maxima, with their prominence, are candidate modes
at each frequency (not yet connected across frequencies). For the
ellipticity angle the resolution is the step of the ellipticity grid
and the kernels are periodic with period 180 degrees (-90 and 90
degrees are the same, purely horizontal motion). The maxima can be
selected by rank and by a phase-velocity interval
(`velocity=(minimum, maximum)`), for instance to follow the
ellipticity of one mode.

## 7. Geometric array response

The response of an array depends only on the station positions. It
tells which wavevectors the array resolves and which it confuses,
independently of the data and of the estimator.

**Array factor and ambiguity.** For stations at $`𝐫_n`$,

```math
B(\Delta𝐤) = \frac{1}{N} \sum_n
    e^{-i\, \Delta𝐤 \cdot 𝐫_n},
\qquad
A(\Delta𝐤) = \left| B(\Delta𝐤) \right|^2
```

$`A(𝟎) = 1`$ and $`A`$ does not change with a translation of
the array. With normalized steering vectors $`𝐚`$,

```math
\left| 𝐚(𝐤_1)^H 𝐚(𝐤_2) \right|^2
    = A(𝐤_2 - 𝐤_1)
```

so $`A`$ is the similarity of two wavevectors as a function of their
difference. Its central region describes the **resolution**; remote
maxima, or connected regions of high response, are **ambiguities**
(aliases, grating lobes of regular arrays, strong sidelobes of
irregular ones). For an irregular array both depend on the direction,
so the full two-dimensional response is the primary result; the
classical rules $`2\pi / d_{\max}`$ (Rayleigh), $`2\pi / d_{\min}`$
(alias period) and $`\pi / d_{\min}`$ (Nyquist) are diagnostic
proxies only.

**Main lobe.** For a level $`T_R`$ (default 0.5) the main lobe is the
connected region of $`A \geq T_R`$ containing the origin; its radius
$`r_R(\theta)`$ is the first crossing of the level along each axis
azimuth, $`2 r_R`$ the full width. The conservative scalar is the
largest radius, the resolution radius $`k_R`$ (`resolution_radius`):
below it two waves are not resolved (a single wave can still be
measured at smaller $`k`$, with a precision limited by noise and
averaging; Section 8). Near the origin
$`A(𝐤) \approx 1 - 𝐤^T 𝐂_r 𝐤`$,
with $`𝐂_r`$ the covariance of the centred coordinates: its
principal directions and the anisotropy ratio describe the geometry.

**Ambiguity.** For a level $`T_A`$ (default 0.5) the regions of
$`A \geq T_A`$ other than the central one are the remote ambiguity
regions; their local maxima are the ambiguity peaks (displacements,
not propagation directions). Two distances must not be confused: the
**peak distance**, from the origin to the nearest significant peak,
and the **onset distance** $`d_{\mathrm{onset}}`$
(`onset_distance`), from the origin to the nearest point of a remote
region; the onset is shorter, because a region crosses the level
before its peak. The peak sidelobe level and the ambiguity margin
$`1 - \max A_{\mathrm{remote}}`$ (a squared correlation, not a
probability) summarize the strength of the ambiguities.

**Search domains.** A search over a domain $`D`$ of wavevectors is
free of ambiguity at level $`T_A`$ if no difference
$`𝐤_2 - 𝐤_1`$ of two of its wavevectors falls into a
remote region. For a full-azimuth disk $`|𝐤| \leq K`$ the
differences fill the disk $`|\Delta𝐤| \leq 2K`$, whence the
conservative limit, the aliasing radius $`k_A`$
(`aliasing_radius`),

```math
k_A = \frac{d_{\mathrm{onset}}}{2}
```

the rigorous form of the factor one half of the classical aliasing
rules (for a regular 3 x 3 array of 10 m: peak at
$`2\pi / 10 = 0.628`$ rad/m, onset 0.531 rad/m, aliasing radius
0.266 rad/m; the limit tends to $`\pi / d`$ as $`T_A`$ tends to one).

**Velocity limits.** At frequency $`f`$ the wavenumber limits become

```math
c_{\max}(f) = \frac{2 \pi f}{k_R}, \qquad
c_{\min}(f) = \frac{2 \pi f}{k_A}
```

geometric bounds at the chosen levels, not confidence limits. When
the aliasing radius is smaller than the resolution radius no
wavenumber is both resolved and unambiguous (e.g. three stations).

**Coordinate uncertainty.** With independent Gaussian errors of
covariance $`𝐒_n`$ and
$`q_n(𝐤) = \exp(-𝐤^T 𝐒_n 𝐤 / 2)`$,
the expected response is

```math
\mathrm{E}\left[A(𝐤)\right]
    = \frac{\left| \sum_n q_n\, e^{-i\, 𝐤 \cdot 𝐫_n}
            \right|^2
            + \sum_n \left(1 - q_n^2\right)}{N^2}
```

where the second term keeps the auto terms, in which the same error
appears twice and cancels. Uncertain positions damp the coherent
structure increasingly with the wavenumber.

A linear array has a valid response but no closed two-dimensional
main lobe (it does not resolve the wavevector across the line): the
limits are refused for it. The numerical grid of the response is
chosen from the proxies only for sampling, and an excessive grid is an
error, never silently coarsened.

## 8. Statistical performance and operational limits

The geometric response is the point-spread function of Bartlett for a
single coherent wave; it is not that of Capon, MUSIC or RTBF, whose
behaviour depends on the SNR, on the number of windows, on the
loading and on the subspace model. $`A = 0.8`$ is not a probability
of aliasing. The performance of an estimator is therefore measured by
Monte Carlo simulation with the estimators themselves:

```math
𝐱_m = s_m\, 𝐚(𝐤_0) + 𝐧_m,
\qquad s_m \sim \mathcal{CN}(0, P_s),
\qquad 𝐧_m \sim \mathcal{CN}(𝟎, P_n 𝐈)
```

(independent realizations, SNR $`P_s / P_n`$ per station). Each trial
builds a sample CSM from $`L`$ realizations, computes the spectrum on
the search grid and classifies its maximum by the error
$`\Delta𝐤 = \hat{𝐤} - 𝐤_0`$ with the
geometric response:

| Outcome | Meaning |
|---|---|
| `BOUNDARY` | maximum on the radial boundary of the search grid |
| `LOCAL` | $`\Delta𝐤`$ in the main lobe (geometrically local, not "accurate") |
| `ALIAS` | $`\Delta𝐤`$ in a remote ambiguity region (an alias family) |
| `OTHER` | anywhere else (noise-driven maxima) |

The probabilities (with Wilson score intervals), the bias and
covariance of the local errors (radial and transverse) and the exact
velocity statistics (each estimate converted,
$`c = 2\pi f / \hat{k}`$) quantify the performance; a sweep over SNR
and windows tells how much data a survey needs.

The **operational assessment** repeats the experiment over a grid of
frequencies, velocities and azimuths and marks a point usable when it
meets explicit criteria (minimum LOCAL probability, maximum ALIAS,
OTHER and BOUNDARY probabilities, maximum relative velocity RMSE,
grid quantization small relative to the resolution radius); otherwise
the first failed criterion is its limiting factor (search boundary,
grid quantization, aliasing, other failures, low SNR, resolution).
There is no opaque quality score. Recommendations explain the limits
and suggest processing changes; the geometry itself is not optimized.

## 9. Frequency domain decomposition

The frequency domain decomposition (FDD; Brincker et al., 2001)
identifies the resonances of a system from its response alone, without
knowing the excitation: the modes of a structure, or the
two-dimensional resonances of a sediment-filled valley recorded by a
linear array (Ermert et al., 2014). It uses the same cross-spectral
matrices as the beamforming (Section 1), with each component of each
station as a channel, but no wave model: no coordinates enter the
computation, and the shapes are read on the geometry of the array.

**Decomposition.** At each frequency the Hermitian, positive
semi-definite CSM of a block is decomposed into its singular values
$`s_1 \geq s_2 \geq \dots \geq s_N`$ and singular vectors
$`𝐮_1, \dots, 𝐮_N`$ (for such a matrix the singular value
decomposition coincides with the eigendecomposition):

```math
𝐂(f) = 𝐔(f)\, 𝐒(f)\, 𝐔^H(f),
\qquad 𝐒 = \mathrm{diag}(s_1, \dots, s_N)
```

For a lightly damped system excited by broadband noise, close to a
natural frequency the response is dominated by one mode: the CSM has
rank close to one, $`s_1`$ peaks and $`𝐮_1`$ is proportional to the
mode shape. The following singular values show closely spaced modes
and the level of the noise. Away from resonances, or for a wavefield
of travelling waves (as in most f-k analyses), several sources
contribute and the vectors are complex and change from block to
block.

**Statistics over blocks.** Each block gives one decomposition. The
singular values are summarized by their geometric mean and the
standard deviation of their natural logarithm, which give the band
$`\bar{s}\, e^{\pm \sigma_{\ln s}}`$. A singular vector is defined up
to a complex factor of unit modulus: each is rotated by
$`-\frac{1}{2} \arg(𝐮^T 𝐮)`$, the angle that maximizes the norm of its
real part, and given the sign that agrees with the previous blocks.
The shape of a mode is the mean of the real parts, with unit norm and
its largest component positive, and its stability is the mean
modal assurance criterion (MAC; Allemang and Brown, 1982) between the
shape and the vectors of the blocks:

```math
\mathrm{MAC}(𝐚, 𝐛) = \frac{|𝐚^H 𝐛|^2}{(𝐚^H 𝐚)(𝐛^H 𝐛)}
```

which is 1 for proportional vectors and 0 for orthogonal ones.

**Peaks.** The candidate resonances are the local maxima of the
geometric mean of $`s_1`$ over the analysis frequencies whose
prominence, the ratio between the peak and the higher of the minima
that separate it from larger peaks, exceeds a threshold (2 by
default): the small fluctuations of the estimate have prominences
close to 1.

**Options.** The windows have a fixed duration (constant resolution,
as in the classical FDD) or a fixed number of cycles (resolution
proportional to the frequency). The matrix can average, besides the
windows of a block, a few sub-frequencies equally spaced in logarithm
around each frequency (`subfrequencies`, relative `bandwidth`), which
adds realizations at the cost of resolution. With
`normalization="coherence"` each matrix is replaced by its coherence,
$`C_{ij} / \sqrt{C_{ii} C_{jj}}`$, multiplied by the mean power
$`\mathrm{tr}(𝐂) / N`$: the differences of gain between the channels
(uncorrected data, different sensors) no longer weigh on the vectors.

The frequency of a peak is limited by the frequency step of the
analysis and by the resolution of the windows, and the FDD gives no
damping: the enhanced FDD gives both.

**Enhanced FDD.** Close to a mode, $`s_1(f)`$ is the auto spectrum of an
equivalent single-degree-of-freedom system (Brincker, Ventura and
Andersen, 2001). Around the peak, the lines whose singular vector
matches the shape at the peak form the bell of the mode:

```math
\mathrm{MAC}\big(𝐮_1(f_{\mathrm{peak}}),\, 𝐮_k(f)\big) \geq \alpha
```

with $`\alpha`$ = `mac` (0.8 by default) and $`k`$ the best matching
singular vector of each line, whose singular value is taken; the bell
is searched on both sides of the peak up to the first line below
$`\alpha`$, and the other lines are set to zero. The inverse Fourier
transform of the bell is the autocorrelation of the free decay of the mode, normalized
to $`r(0) = 1`$:

```math
r(\tau) \approx e^{-\zeta \omega_n \tau} \cos(\omega_d \tau),
\qquad \omega_d = \omega_n \sqrt{1 - \zeta^2}
```

Its extrema (peaks and troughs, refined by a parabola through three
samples) are half a damped period apart, so that their amplitudes
$`|r_k|`$ decrease as $`e^{-k \delta / 2}`$, with $`\delta`$ the
logarithmic decrement. A linear regression of $`2 \ln |r_k|`$ on
$`k`$, on the extrema whose amplitude lies between the `levels`
(0.1 and 0.9 by default, stopping at the first one below the lower
level, where noise dominates), gives $`\delta`$ as the slope; a
regression of the zero crossings (linear interpolation) on their
index gives half the damped period. Then

```math
\zeta = \frac{\delta}{\sqrt{4 \pi^2 + \delta^2}},
\qquad f_n = \frac{f_d}{\sqrt{1 - \zeta^2}}
```

The inverse transform needs equally spaced lines of a single spectral
estimate: the enhanced FDD computes, around each mode, a local
spectrum with windows of fixed duration $`T`$ (averaging of all the
windows as in Welch, 1967, Hann taper, 50 % overlap) on the FFT lines
$`k / T`$. Windows with a fixed number of cycles, valid for the FDD,
would give each line its own resolution and leakage, and a bell that
is not the spectrum of one correlation function. $`T`$ is a number of
cycles of the mode (`cycles`, 100 by default: about 50 lines inside
the half-power bandwidth $`2 \zeta f_n`$ at 2 % damping) or a given
duration. Leakage widens the bell and overestimates the damping; a
bell truncated by a high $`\alpha`$ underestimates it; the number of
windows sets the statistical scatter.

# Part II. Use

```python
from shakelab.exploration.yarray import ArrayData, Beamformer
```

The package exposes the main classes and functions; the modules are
`data`, `config`, `estimators`, `beamforming`, `fdd`, `statistics`,
`response`, `performance`, `operational`, `io`, `monitor` and `cli`,
the figures `shakelab.graphics.yarray`. Options with fixed values are
enumerations whose string values are accepted everywhere
(`method="capon"`).

## 10. Data

An `ArrayData` holds the stations (`StationData`: records, ENU
coordinates, orientation of the components, optional coordinate
covariance and geographic position). It is filled in one go or
progressively, from any source:

| Source | Use |
|---|---|
| Array file (JSON: local coordinates, one waveform file per component, reader options) | `ArrayData("array.json")`, or `array.read("array.json")` |
| Waveforms read by ShakeLab (StreamCollection, Stream, records: any format and source), with metadata or a table of coordinates | `ArrayData(collection, select=..., coordinates=..., coordinate_covariance=...)`, or `array.add(...)` |
| Stations already built | `array.add(station)` |

```python
from shakelab.exploration.yarray import ArrayData

array = ArrayData()
array.read("array.json")
array.add(records, coordinates={"XX.A10.": (120.0, -35.0)})
array.remove("XX.A03.")
array.info()
```

The stations need not share their time axis: each keeps the interval
covered by all its channels (or `starttime`, `endtime`), and the
window common to all of them (`starttime`, `endtime`, `nsamp` of the
array) is taken when the array is processed (`array.window(...)`,
used by the beamformer). Stations with gaps are excluded with an
`ArrayDataWarning` (`gaps="drop"`) or raise (`gaps="raise"`); samples
that are not simultaneous are an error (no implicit resampling);
sampling interval and units must be the same for all the stations.

With geographic positions in the metadata the coordinates are
computed in a local ENU frame (`ArrayData.frame`). Its origin is fixed:
the `reference` given to `ArrayData` or written in the array file
(`"reference": [lon, lat, elevation]`), otherwise the position of the
first station added; adding stations does not move it.
`set_reference(point)` (or `set_reference()`, the mean of the current
stations) changes it and recomputes the coordinates.

Positions are optional. A station without a metadata location and
without explicit coordinates (or, in the array file, without
`coordinates`) is added without a position (`station.has_position`,
`array.has_positions`). The wavenumber methods need every position:
`array.coordinates`, the `Beamformer`, `ArrayGeometry` and
`ArrayResponse` raise a `ValueError` naming the stations without one.
The frequency domain decomposition does not use them; its results keep
the positions that are known (NaN for the others), to draw the mode
shapes on the geometry of the array.

### Array file

An array file describes an array whose recordings have no usable
metadata (for example SAC files of a temporary survey): the stations,
their local coordinates and the waveform file of each component. It is
a JSON object (version 1):

```json
{
  "version": 1,
  "coordinate_system": {"type": "local_cartesian", "unit": "m"},
  "reference": [13.25, 46.06, 120.0],
  "data_root": "sac",
  "reader": {"format": "sac", "options": {}},
  "stations": [
    {"id": "S1001",
     "coordinates": {"x": 2048.0, "y": 1964.0, "z": 0.0},
     "coordinate_covariance": [[0.01, 0.0], [0.0, 0.01]],
     "components": {"Z": "S1001.Z.sac", "N": "S1001.N.sac",
                    "H1": {"file": "S1001.1.sac",
                           "orientation": [0.7, 0.7, 0.0],
                           "reader": {"format": "sac"}}}}
  ]
}
```

| Field | Content |
|---|---|
| `version` | format version, 1 (required) |
| `coordinate_system` | `local_cartesian`, unit `m` |
| `reference` | optional geographic position (longitude, latitude[, elevation in m]) of the origin of the coordinates |
| `data_root` | folder of the waveform files (default: the folder of the array file) |
| `reader` | optional arguments of `shakelab.signals.io.reader` (`format`, `byte_order`, `options`) for all the files |
| `stations` | list of stations (required) |
| `id` | station name |
| `coordinates` | East, North[, Up] in metres, as an object `x`, `y`, `z` or a list |
| `coordinate_covariance` | optional covariance of the horizontal coordinates, in square metres |
| `components` | waveform file of each component: a file name, or an object with `file`, `orientation` and `reader` |

Components named Z, N and E (also UP, VERTICAL, NORTH, EAST, X, Y, and
1, 2, 3 for E, N, Z as in the ArrayTool3C file lists) imply their
orientation; other names need an explicit `orientation` (direction of
the sensor axis as east, north, up components). Each file must contain
one record; the components of a station are cut to their common time
window. An array file takes no reader options in `ArrayData` (they are
written in the file).

## 11. Configuration

`BeamformingConfig` collects every parameter. Each one can be given
singly, with the name of the corresponding command-line option
(`--grid-max` is `grid_max`), or grouped in the objects
`FrequencyConfig`, `WindowConfig`, `SearchGrid` and `EllipticityGrid`
(immutable, also as dictionaries with the fields of the JSON
sections); single parameters apply on top of the objects and of the
file. The configuration changes with `update` (all the values are
validated together) or by assigning a parameter.

```python
from shakelab.exploration.yarray import BeamformingConfig

config = BeamformingConfig(fmin=2.0, fmax=12.0, nfreq=30,
                           grid_min=100.0, grid_max=2000.0, grid_size=100,
                           window_type="fixed", window_length=5.0)
config = BeamformingConfig("capon.json", method="music")  # file + changes
config.update(components="radial", maxima=3)
config.window_length = 4.0
config.write("radial.json")
```

| Parameter | Default | Meaning |
|---|---|---|
| `frequencies` | required | `FrequencyConfig`, dictionary or explicit values; or `fmin`, `fmax`, `nfreq`, `frequency_scale` |
| `grid` | required | `SearchGrid` or dictionary; or `grid_min`, `grid_max`, `grid_size`, `grid_quantity`, `grid_scale`, `azimuth_size`, `azimuths` |
| `window` | `WindowConfig()` | windows and covariance blocks; or `window_type`, `window_length`, `window_taper`, `window_overlap`, `windows_per_block`, `block_overlap` |
| `method` | `capon` | `bartlett`, `capon`, `music` |
| `components` | `vertical` | `vertical`, `radial`, `transverse`, `three_component_separate`, `rayleigh_joint`, `rayleigh_angular`, `rayleigh_direct` |
| `diagonal_loading` | 0.001 | relative loading of Capon |
| `music_sources` | 1 | sources excluded from the MUSIC noise subspace |
| `ellipticity_grid` | -85 to 85 degrees, 61 | `EllipticityGrid` of RTBF; or `ellipticity_angle_min`, `ellipticity_angle_max`, `ellipticity_angle_size` |
| `spectra` | False | keep also the complete spectra of every block (the maxima are always kept); they can take much memory |
| `maxima` | 1 | maxima kept per spectrum |
| `peak_neighborhood` | (3, 3) | neighbourhood of a maximum (azimuth, radial) |
| `minimum_relative_power` | 0 | smallest power of a maximum, relative to the largest |
| `statistical_uncertainty` | `none` | `empirical` (Bartlett, Capon; single components) |
| `geometric_uncertainty` | `none` | `expected_coherence` (Bartlett; needs coordinate covariances) |
| `strict` | True | more windows than the matrix dimension for Capon and MUSIC |
| `transform` | `direct` | `direct`, `rfft` (fixed windows) |
| `fft_selection` | `nearest` | `nearest`, `interp` (with `rfft`) |
| `fft_padding` | None | zero padding (1 for nearest, 8 for interp) |

| `FrequencyConfig` | Default | Meaning |
|---|---|---|
| `values` | None | explicit frequencies (Hz) |
| `minimum`, `maximum`, `size` | None | range of `size` frequencies; with the RFFT a range without `size` takes every bin |
| `scale` | `log` | `log`, `linear` |

| `WindowConfig` | Default | Meaning |
|---|---|---|
| `type` | `cycles` | `cycles` (length in cycles of each frequency) or `fixed` (seconds) |
| `length` | 20 | cycles or seconds |
| `taper` | 0.1 | Tukey fraction |
| `overlap` | 0.5 | overlap of the windows |
| `windows_per_block` | 40 | windows of a covariance block |
| `block_overlap` | 0.5 | overlap of the blocks |

| `SearchGrid` | Default | Meaning |
|---|---|---|
| `minimum`, `maximum`, `size` | required | radial range and samples |
| `quantity` | `velocity` | `velocity`, `slowness`, `wavenumber` |
| `scale` | log (linear for wavenumbers) | `log`, `linear` |
| `azimuth_size` | 180 | uniform azimuths over 360 degrees |
| `azimuths` | None | explicit azimuths |

The JSON file has the sections `frequencies`, `spectral` (`transform`,
`fft_selection`, `fft_padding`), `method` (`name`, `diagonal_loading`,
`music_sources`), `components`, `grid`, `window`, `ellipticity`
(`minimum_angle`, `maximum_angle`, `size`), `output` (`spectra`,
`maxima`, `neighborhood`, `minimum_relative_power`; an older `mode`
key is read as `spectra`), `uncertainty` (`statistical`,
`geometric`) and `strict`:

```json
{
  "version": 1,
  "frequencies": {"minimum": 2.0, "maximum": 12.0, "size": 100,
                  "scale": "log"},
  "spectral": {"transform": "rfft", "fft_selection": "nearest"},
  "method": {"name": "capon", "diagonal_loading": 0.001},
  "components": "vertical",
  "grid": {"quantity": "velocity", "scale": "log", "minimum": 100.0,
           "maximum": 2000.0, "size": 100, "azimuth_size": 180},
  "window": {"type": "fixed", "length": 5.0, "taper": 0.1,
             "overlap": 0.5, "windows_per_block": 40,
             "block_overlap": 0.5},
  "output": {"spectra": true, "maxima": 3},
  "uncertainty": {"statistical": "none", "geometric": "none"},
  "graphics": {"dispersion": true, "ranks": [1, 2, 3]}
}
```

Missing fields take the defaults; unknown sections are ignored by the
processing, so the section `graphics` of the command line can live in
the same file.

## 12. Processing and results

```python
from shakelab.exploration.yarray import ArrayData, Beamformer

array = ArrayData("array.json")
beamformer = Beamformer(array, "capon.json")
result = beamformer.process()                       # common window
part = beamformer.process("2026-11-01T10:00", "2026-11-01T10:10")
result.info()
result.dispersion().info()

for maximum in result.maxima_at(rank=1):
    print(maximum.time, maximum.frequency, maximum.velocity,
          maximum.azimuth, maximum.power)
```

`Beamformer(array, config)` keeps the array by reference (stations
added later are processed too) and a copy of the configuration
(`beamformer.config`, which can be changed between runs).
`cache=True` keeps the steering vectors while station coordinates and
grid do not change (repeated analyses, real time). Frequencies whose
windows do not fill one block give no result and a
`BeamformingWarning`, not an error.

`BeamformingResult` holds the maxima (`SpectralMaximum`), the spectra
if requested (`spectrum(frequency, component)` returns a
`FrequencySpectrum` with the power of every block on the grid, and
its standard error or ellipticity), the configuration, the stations
(`sids`, `coordinates`, and `stations_at(frequency, block)` for the
stations of each block), the block times and the units of the power
(`power_units`). `dispersion(rank=1)` gives, at each frequency, the
median velocity of the maxima with its 16-84 % interval, their mean
azimuth and their number (`DispersionCurve`). `append` adds the
results of later windows (blocks renumbered; results with the same
configuration, frequencies and units).

```python
from shakelab.exploration.yarray import BeamformingResult

result.write("daily.jsonl")
beamformer.process(start, end).write("daily.jsonl", append=True)
daily = BeamformingResult("daily.jsonl")              # joined
```

A result file holds one result, possibly written in several pieces:
`write(append=True)` checks that the configuration is the same (see
"Result files" below for the format).

The density of the maxima:

```python
density = result.density(coordinate="log_velocity",
                          weighting="relative_power")
summary = density.summary()              # mode, median, quantiles
peaks = density.peaks(3, minimum_relative_density=0.1)
ellipticity = result.density(coordinate="ellipticity_angle",
                             weighting="relative_power",
                             velocity=(150, 400))     # one mode (Rayleigh)
```

Figures (`shakelab.graphics.yarray`, matplotlib imported on use):
`plot_dispersion`, `plot_ellipticity`, `plot_spectrum` (power,
standard error or ellipticity; wavenumber plane or azimuth and grid
coordinates), `plot_density`. See the examples `01_array.py` to
`10_density.py` in
[`examples/exploration/yarray/`](../../examples/README.md#exploration).

### Frequency domain decomposition

```python
from shakelab.exploration.yarray import ArrayData, FDD, FDDConfig

array = ArrayData("array.json")
config = FDDConfig(fmin=0.1, fmax=1.0, nfreq=100,
                   frequency_scale="linear", window_length=20,
                   windows_per_block=20, components="transverse",
                   azimuth=30.0, normalization="coherence")
result = FDD(array, config).process()
result.info()

frequency, mean, log_std = result.singular_values()
peaks = result.peaks(fmin=0.2, fmax=0.8)
mode = result.mode(peaks[0])
mode.info()
```

`FDDConfig` has the frequencies and windows of `BeamformingConfig`
(same names; with fixed windows and no `nfreq`, the FFT frequencies of
the window between `fmin` and `fmax`), plus `subfrequencies`,
`bandwidth`, `normalization` (`"none"` or `"coherence"`),
`components` (`"vertical"`, `"north"`, `"east"`, `"radial"` and
`"transverse"` with `azimuth`, `"horizontal"`, `"three_component"`)
and `vectors`, the number of singular vectors kept for each block (3).
`FDD(array, config)` keeps the array by reference and a copy of the
configuration, like the `Beamformer`; frequencies without complete
blocks give an `FDDWarning`.

`FDDResult` holds, by frequency, the block times, all the singular
values and the first singular vectors of every block (with normalized
phase and sign), the channels (`stations`, `components`, `channels`),
their coordinates and the units (`power_units`).
`singular_values()` returns the geometric mean and the logarithmic
standard deviation; `peaks(fmin, fmax, index, prominence, count)` the
frequencies of the peaks; `mode(frequency, index)` a `Mode` with the
shape, its standard deviation among the blocks, the singular value,
the mean MAC and the channels. `mac(a, b)` compares two shapes. Files
work as for the beamforming: `write(..., append=True)` accumulates
results with the same configuration and channels, read back joined by
`FDDResult(file)`.

The enhanced FDD (Section 9) is a method of the processor, which has
the data: it computes its own local spectrum around each mode.

```python
fdd = FDD(array, config)
result = fdd.process()
modes = fdd.efdd(result.peaks(), mac=0.8)  # cycles=100 or duration=...
for mode in modes:
    mode.info()                            # frequency, damping, shape
```

`efdd(frequencies, mac, cycles, duration, levels, starttime, endtime,
gaps)` returns one `Mode` for each frequency (the peak is searched
within 5 %), with `method="efdd"`, the natural frequency, the damping
ratio, the shape at the peak and the mean MAC of the bell; `mode.efdd`
(`EFDDDetails`) keeps the local spectrum, the bell, the correlation,
its extrema and crossings and the regression, for checking the
estimate. The components and azimuth of the configuration are used,
not its windows or normalization (the bell must be a power spectrum).

Figures: `plot_singular_values(result, count=3)` (geometric means with
their bands and the peaks), `plot_mode_shape(mode, along=...)` (shape
along the main direction of the stations, a coordinate, `"z"` with
the height on the vertical axis for a building, or the station index
when positions are missing), `plot_mac(modes)` (MAC between modes)
and `plot_efdd(mode)` (bell and fitted decay). See
[`17_fdd.py`](../../examples/exploration/yarray/17_fdd.py) and
[`19_efdd.py`](../../examples/exploration/yarray/19_efdd.py).

### Result files

Results are written as JSON Lines (`.jsonl`): text, one JSON object per
line, readable by any language. The first line is a header with the
kind of result (`"beamforming"` or `"fdd"`), the complete
configuration, the frequencies and the units (for the FDD also the
channels and their coordinates); each result written adds one line.
Appending (`write(file, append=True)`, the `output` of the monitors,
`slab-yarray --append`) only adds a line, also over days of real-time
processing, and is refused if the header differs. The name `.json` is
refused, because the file is not one JSON document.

| Result | In the `.jsonl` line | In the `.npz` file (`arrays=True`) |
|---|---|---|
| Beamforming | block times, stations and coordinates of each block, table of the maxima (block, time, frequency, component, rank, wavenumber, azimuth, velocity, slowness, power, standard error, ellipticity and its angle) | the spectra on the grid, if the configuration keeps them (`spectra=True`) |
| FDD | block times, first singular value and first singular vector (real and imaginary parts) of every block | all the singular values and the kept singular vectors (`vectors` of the configuration) of every block |

The arrays of each line go to a compressed `.npz` file in the folder
`<name>_arrays` next to the result file, named after the first block
time (`fdd_arrays/20261116T101500.npz`) and read without pickle. They
are written by default (`write(..., arrays=False)` writes the line
only) and read when every line of the file has them;
`FDDResult(file, arrays=False)` reads only the lines. The line alone is
enough for the usual work: the maxima, the dispersion curves and the
density of the beamforming; the singular-value spectrum, its peaks and
the first mode of the FDD. A result read without arrays keeps only the
first singular value and vector, and is not joined with a complete
one.

```python
result.write("fdd.jsonl")                    # line + fdd_arrays/*.npz
result.write("fdd.jsonl", append=True)       # one more line
light = FDDResult("fdd.jsonl", arrays=False)
result.to_csv("fdd.csv")                     # optional table
```

`to_csv` writes an optional table for spreadsheets: the maxima of the
beamforming, or the first singular value and vector of every block of
the FDD (missing values empty, times in ISO 8601). The result file is
the complete record.

## 13. Command line

`slab-yarray` (installed with ShakeLab; `python -m
shakelab.exploration.yarray` is equivalent) runs the beamforming and
the array response; the frequency domain decomposition is used from
Python and in real time.

```sh
slab-yarray --input array.json --config capon.json -o result.jsonl \
    --plot-dispersion --save-figures figures
slab-yarray --input array.json --config capon.json --method music \
    --fmin 2 --fmax 8 --nfreq 20 -o music.jsonl     # options replace the JSON
slab-yarray --data "data/*.mseed" --metadata array.xml \
    --select "XX.*..HH?" --config rtbf.json -o rtbf.jsonl --append
slab-yarray --result rtbf.jsonl --plot-spectrum --frequency 4.5 --show
slab-yarray response --input array.json --fmin 1 --fmax 20 --plot-all
```

| Group | Options |
|---|---|
| Input (one) | `--input` array file; `--data` waveform files with `--metadata`, `--select`, `--start`, `--end`; `--result` saved result to plot |
| Configuration and output | `--config`, `--write-config`, `-o/--output` (result file, `.jsonl`), `--append`, `--no-arrays` (no `.npz` spectra), `--csv` (table of the maxima, also of a `--result`), `--owrite` |
| Processing (replace the JSON) | `--frequencies` or `--fmin --fmax --nfreq`, `--frequency-scale`, `--method`, `--components`, `--diagonal-loading`, `--music-sources`, `--grid-quantity`, `--grid-scale`, `--grid-min`, `--grid-max`, `--grid-size`, `--azimuth-size`, `--window-type`, `--window-length`, `--window-taper`, `--window-overlap`, `--windows-per-block`, `--block-overlap`, `--ellipticity-angle-min/max/size`, `--spectra`, `--maxima`, `--minimum-relative-power`, `--statistical-uncertainty`, `--geometric-uncertainty`, `--strict`, `--transform`, `--fft-selection`, `--fft-padding` |
| Figures | `--plot-dispersion`, `--plot-ellipticity`, `--plot-spectrum` (`--frequency`, `--block`, `--spectrum-view`, `--spectrum-quantity`), `--plot-density` (`--density-coordinate`), `--component`, `--ranks`, `--show`, `--save-figures`, `--figure-format`, `--dpi` |

The figures can also be requested in the section `graphics` of the
configuration file (`dispersion`, `ellipticity`, `component`, `ranks`,
`spectrum` {`enabled`, `frequency`, `block`, `view`, `quantity`},
`density` {`enabled`, `coordinate`}, `show`, `save` {`directory`,
`format`, `dpi`}); the options of the command line prevail.

`slab-yarray response` analyses the geometry only (Section 7): input
`--input`, `--data` or `--coordinates` (a text file of lines `x y [z]`
or `name x y [z]`, to check a layout before deployment); levels
`--resolution-level`, `--ambiguity-level`; `--expected` (coordinate
uncertainty of the array); grid `--grid-max --grid-size`; velocity
limits at `--frequencies` or `--fmin --fmax --nfreq`; `--search-max K`
to check a search domain; figures `--plot-geometry`,
`--plot-response`, `--plot-resolution`, `--plot-velocity`,
`--plot-search-domain` or `--plot-all`. The Monte Carlo analysis
needs too many choices for a command line and is used from Python.
See
[`11_command_line.sh`](../../examples/exploration/yarray/11_command_line.sh).

## 14. Real time

```python
from shakelab.exploration.yarray import BeamformingMonitor
from shakelab.signals.io import metadata_reader
from shakelab.signals.seedlink import Client

client = Client("server", metadata=metadata_reader("array.xml"))
client.select("XX.A*..HH?")
client.set_buffer(max_duration=900)
monitor = BeamformingMonitor(config, window=600, output="monitor.jsonl")
monitor.on_result(lambda monitor, result: print(result.summary()))
client.every(600, monitor)        # clock="data" to replay an archive
client.start()
...                               # from any thread, while it runs:
monitor.status(); monitor.latest(); monitor.result(start, end)
client.stop(); monitor.close()
```

The monitor is called by the client in its receiving thread: it only
copies the latest `window` seconds of the buffer (`client.get_data`)
and queues them for its worker thread, so the reception is never
delayed. The queue is short (`queue_size=1`): when the processing is
slower than the calls the oldest pending window is dropped and
reported (`on_skip`), so the monitor never accumulates delay. Windows
that cannot be processed (buffer shorter than the window, too few
stations, no complete block) are skipped and reported; an unexpected
error stops the worker and is raised by `close()`. `status()`,
`latest()`, `results(start, end)`, `result(start, end)` (joined) and
`array()` can be called at any time and return copies; each result can
also be appended to a file (`output`, with the arrays unless
`arrays=False`), and every analysis, skipped window or error to a log
(`log`, JSON Lines: wall time, event, window, duration and stations,
or reason).
`background=False` processes in the call;
`monitor.process(collection, start, end)` processes given data. Any
object with `get_data(select)` returning a `StreamCollection` can feed
the monitor. See
[`12_realtime.py`](../../examples/exploration/yarray/12_realtime.py).

There is one monitor for each processing, with the same options and
queries: `BeamformingMonitor` (a `Beamformer` with the cache of the
steering vectors; results `BeamformingResult`) and `FDDMonitor` (an
`FDD`; results `FDDResult`). Both derive from `ArrayMonitor`, which
holds the common machinery and is not used directly.

```python
from shakelab.exploration.yarray import FDDMonitor

monitor = FDDMonitor(config, window=1800, output="fdd.jsonl")
client.every(300, monitor)
...
monitor.stations                  # fixed by the first result
monitor.result().mode(1.2).info()
```

Results of the decomposition join only with the same channels, so an
`FDDMonitor` always processes the same stations: those given with
`stations=[...]` (in that order) or else those of its first result.
Other stations of the buffer are ignored; a window in which one of
them is missing is skipped. By default (`gaps="raise"`) a window in
which one of them has gaps is skipped too, rather than processed
without it. The stations need no position. See
[`18_realtime_fdd.py`](../../examples/exploration/yarray/18_realtime_fdd.py).

## 15. Survey planning

```python
from shakelab.exploration.yarray import (
    ArrayData, ArrayResponse, PerformanceConfig,
)

response = ArrayResponse(ArrayData("array.json"))  # or coordinates
analysis = response.analyze()                      # levels 0.5
analysis.info()
limits = analysis.velocity_limits([1, 2, 5, 10, 20])
response.search_domain(0.1).info()

config = PerformanceConfig(snr_db=0.0, method="capon", grid_min=0.0,
                           grid_max=0.2, grid_size=41, trials=500,
                           realizations=40, seed=1)
wave = {"velocity": 400.0, "frequency": 8.0, "azimuth": 30.0}
response.performance(config, **wave).info()
response.performance_sweep(config, snr_db=[-10, 0, 10],
                           realizations=[10, 20, 40], **wave).info()
assessment = response.assess(config, frequencies=[2, 4, 8, 16],
                             velocities=[150, 300, 600, 1200])
assessment.info()
```

`ArrayResponse` accepts an `ArrayData` (and then its coordinate
covariances, for the expected response: `analyze(expected=True)`), an
`ArrayGeometry` or coordinates, and is computed once: the response does
not follow later changes of the array. The wavenumber grid is chosen
from the geometry (extent 1.25 times $`2\pi / d_{\min}`$, 12 samples
per $`2\pi / D`$, with $`d_{\min}`$ the minimum spacing and $`D`$ the
maximum aperture); `grid_max` and `grid_size` replace either choice.

`PerformanceConfig` works as `BeamformingConfig` (single parameters,
JSON file in the constructor, `update`, `copy(**changes)`, `read`,
`write`): `snr_db` and the search grid (`grid_min`, `grid_max`,
`grid_size`, `azimuth_size` or `azimuths`; always of wavenumbers, the
experiment does not depend on the frequency) are required; `method`
(Capon by default), `diagonal_loading`, `music_sources`, `strict`;
`trials`, `realizations`, `seed` (reproducible runs; point `i` of a
sweep or assessment uses a derived seed); `resolution_level`,
`ambiguity_level` and `boundary_margin` of the classification. It is
not derived from a beamforming configuration, whose grid is usually of
velocities and whose overlapping windows are fewer independent
realizations than `windows_per_block`.

The plane wave is given by `wavenumber`, or by `velocity` and
`frequency`, with `azimuth` (0) and `signal_power` (1); a known
frequency is kept in `result.source` and used by
`result.velocity_statistics()`. The criteria of `assess` are
parameters: `min_local` (0.9), `max_alias` (0.05), `max_other`
(0.05), `max_boundary` (0.01), `max_velocity_rmse` (0.1, relative) and
`max_grid_error` (0.1, relative to the resolution radius), with
`conservative_azimuth` (all the azimuths, or at least half). Figures:
`plot_array_geometry`, `plot_array_response`, `plot_resolution`,
`plot_search_domain`, `plot_velocity_limits`, `plot_performance`,
`plot_performance_sweep`, `plot_operational`. See
[`13_response.py`](../../examples/exploration/yarray/13_response.py),
[`14_performance.py`](../../examples/exploration/yarray/14_performance.py)
and [`15_operational.py`](../../examples/exploration/yarray/15_operational.py).

## 16. Practical guidance and processing times

1. Check the array first (`slab-yarray response`): the velocity limits
   tell which frequencies and velocities the grid should cover.
2. Keep the configuration in a version-controlled JSON file and start
   with a compact run (`spectra=True`, few frequencies) to inspect
   spectra, number of blocks and windows per block.
3. For Capon and MUSIC use more windows per block than stations
   (`2N` for RTBF); prefer more windows to more loading.
4. Use the direct transform for cycle windows or arbitrary
   frequencies; the RFFT with nearest bins for many frequencies on
   fixed windows (the frequencies become the bins); the interpolated
   RFFT keeps the frequencies but is not faster than the direct
   transform when beamforming dominates.
5. The Rayleigh processing works on matrices of dimension `2N`
   (`3N` for ARDS) and needs more windows per block; the RTBF also
   scans an ellipticity grid, while ARTBF and ARDS find the angle
   exactly and are the cheapest of the three.
6. For the FDD use fixed windows long enough to separate the modes,
   and the coherence normalization when the channels have different
   gains (uncorrected data, different sensors). Before reading a peak
   of the first singular value as a mode, check its shape and its
   mean MAC: close to 1 for a standing resonance, low for travelling
   waves.
7. Use the enhanced FDD on well separated modes and long records (its
   damping scatters by about 10 % with 20 minutes of data), and check
   each estimate with `plot_efdd`.
8. Save the results and plot or post-process them from the file,
   rather than repeating the processing.

Times on the nine SESAME stations (three components, 60 s;
[`16_profiling.py`](../../examples/exploration/yarray/16_profiling.py),
one core of a laptop, indicative only): vertical Capon on
180 azimuths x 100 radial values, fixed 2 s windows, takes about
0.5 s for 50 frequencies with the direct transform and 0.13 s with the
nearest-bin RFFT (which keeps 15 distinct bins), and grows linearly
with the frequencies for the direct and the interpolated transforms;
for 10 frequencies, RTBF with 61 ellipticities takes about 0.4 s, ARTBF
and ARDS about 0.25 s; the cache of the steering vectors halves the
time of repeated analyses of the same array. With the 38 stations of
the full simulation, single-component processing takes of the order of
15 s (Bartlett, 80 spectra), the Rayleigh processing about 1.4 s for
RTBF, 1.2 s for ARTBF and 0.9 s for ARDS at 4 frequencies with 130
windows per block.

The decomposition is much cheaper than the beamforming, since it has no
search grid: on the nine SESAME stations, 50 frequencies with 10
windows per block take about 0.02 s for the vertical component and
0.04 s for the three components (27 channels) with fixed 5 s windows,
and about 2.5 times more with windows of 20 cycles. On the four-storey
frame of Section 18 (four channels, 20 minutes at 50 Hz, 20 s windows,
191 frequencies) the FDD takes about 0.14 s and the enhanced FDD of
the four modes 0.13 s with windows of 100 cycles (0.3 s with 200).

# Part III. Design, validation and limits

## 17. Design choices

- **Layers.** Generic primitives in `shakelab.signals` (no array
  concepts) < `data` < `estimators` <
  `config` < processing (`beamforming` and `statistics`, `fdd`; for
  planning `response`, `performance` and `operational`) < applications
  (`monitor`, command line, graphics). The numerical modules do no
  plotting and keep no global state; files are read and written by
  `io` (array and result files), called by the configuration and
  result classes.
- **Array data filled progressively.** The stations keep their own
  time axis and the common window is taken when processing
  (`ArrayData.window`, used by `Beamformer`), so that stations from
  different sources can be added one at a time. `StationData` is
  immutable and the builders create new records: later changes of the
  input (for example a SeedLink buffer) do not affect the array.
- **Fixed origin.** The origin of the ENU frame is the first station
  added with a geographic position, or the given `reference`; adding
  stations does not move it, `set_reference` recomputes the
  coordinates.
- **No implicit corrections.** Stations not covering the window
  without gaps are excluded with a warning (`gaps="drop"`) or raise an
  error (`gaps="raise"`); samples not simultaneous between stations are
  an error (no implicit resampling).
- **One name per parameter.** The parameters of `BeamformingConfig`
  have the names of the command-line options; the grouped objects
  (`WindowConfig`...) are immutable and the configuration changes with
  `update`, validated as a whole. `Beamformer`, `FDD`, the monitors
  and the results keep their own copy. Options with fixed values are
  enumerations, whose string values are accepted everywhere.
- **Results for real time.** A result file holds one result, possibly
  written in pieces (`append=True`, with the same configuration): JSON
  Lines, so that appending adds one line without rewriting the file;
  the configuration is in the header, so that every file is
  self-contained. The bulky arrays (spectra, all the singular values
  and vectors) go to separate `.npz` files, optional, so that the text
  file stays light. The stations and coordinates are recorded for each
  block, so that results of windows with different stations join
  correctly. Frequencies whose windows do not fill a covariance block
  give no result and a warning instead of an error, as needed for
  short windows; `Beamformer(cache=True)` keeps the steering vectors
  while coordinates and grid do not change. The monitors copy the
  buffer in the receiving thread and process in their own, with a
  short queue that drops the oldest window rather than accumulating
  delay.
- **Planning from a snapshot.** `ArrayResponse` is computed once from
  the coordinates of the moment (an `ArrayData` changed later does not
  change it) and is the entry point of every planning analysis.
  `PerformanceConfig` shares the mechanics of `BeamformingConfig` but
  is not derived from it: the Monte Carlo grid is of wavenumbers,
  independent of the frequency, and overlapping windows are fewer
  independent realizations than `windows_per_block`.
- **FDD on the shared matrices.** The cross-spectral matrices by
  frequency and block are a generic object of `shakelab.signals`
  (`crossspectral.CrossSpectra`), with windows in cycles (blocks that
  change with the frequency) and sub-frequency averaging; the FDD
  adds the decomposition, the alignment of the vectors and the
  statistics. `FDD` and `FDDResult` follow `Beamformer` and
  `BeamformingResult` (array by reference, copy of the configuration,
  files with `append`); FDD results join only with the same channels,
  since their statistics combine the blocks channel by channel. The
  result files hold beamforming or FDD results, each class reading
  only its own.
- **Enhanced FDD on its own spectrum.** The bell must be transformed
  to a correlation, which needs equally spaced lines of one estimate:
  the enhanced FDD is a method of `FDD`, which has the data, and
  computes a local spectrum with fixed windows for each mode (a number
  of cycles of the mode, like the windows in cycles, frozen at the
  peak), whatever the windows of the FDD. Its details stay with the
  mode (`Mode.efdd`), so that every estimate can be checked
  (`plot_efdd`).

## 18. Validation

YARRAY derives from development code written in August 2026 for the
MODERN project, brought into ShakeLab with tests against the original
results:

- On the SESAME simulation (nine three-component stations, 60 s) the
  maxima (wavenumber, azimuth, rank, block, ellipticity) are identical
  to those of the August code for every estimator and component mode
  (`tests/exploration/yarray/data/reference.json`).
- Densities, summaries and peaks of the maxima agree within 1e-14,
  except for one correction: the peak prominence is measured from the
  saddle towards a higher peak (`scipy.signal.find_peaks`), not from
  the lowest density on each side, which overrated secondary modes on
  the flank of a stronger one.
- The planning results (responses, limits, Monte Carlo draws and
  outcomes, envelopes) are identical bit for bit with the same seeds.
- The sign of the ellipticity is checked on the fundamental Rayleigh
  mode of a homogeneous half-space (retrograde, $`e = -0.68`$ for a
  Poisson solid, with the theoretical $`c / \beta`$ and H/V ratio), and
  the closed-form ellipticity angle of ARTBF and ARDS against a fine
  grid search.

The frequency domain decomposition is checked on a four-storey shear
frame excited by white noise at every floor, simulated exactly by modal
superposition (2 % damping, natural frequencies from 1.5 to 8.1 Hz, 20
minutes): with 20 s windows (FFT lines every 0.05 Hz) the peaks of the
first singular value fall within 1.5 lines of the natural frequencies,
the shapes have MAC above 0.999 with the exact modes and the first
singular value exceeds the second by more than 20 times at every mode,
with several noise realizations. The cross-spectral matrices agree with
the theoretical matrix of mixed white noise, and the coherence
normalization makes the vectors independent of the channel gains.

The enhanced FDD is checked on the same frame (2 % damping, 20
minutes, windows of 100 cycles): over six noise realizations the
natural frequencies are within 0.3 % on average (scatter 0.1-0.3 %)
and the damping ratios between 2.0 % and 2.2 % on average (scatter
0.1-0.3 % of damping); windows of 200 cycles give 1.8-2.0 % with a
larger scatter (fewer windows), and the MAC threshold (0.8 or 0.9)
changes the estimates little on well separated modes. The test
requires 0.5 % in frequency and 0.6 % of damping on one realization.

Differences from the August code that change numbers: the power is in
physical units (power spectral density), $`N \Delta / 2`$ times the old
value (squared for the RTBF objective), with maxima and relative
thresholds unchanged; the defaults of the windows are those of the
August JSON files (overlap 0.5, 40 windows per block, block overlap
0.5), not those of the August Python classes (0, 20, 0); the default
estimator of the planning is Capon (Bartlett before).

## 19. Current limitations and planned extensions

Current limitations:

- Capon estimated from $`L`$ windows underestimates the power by about
  $`(L - n + 1) / L`$ for a matrix of dimension $`n`$ (Capon and Goodman,
  1970): -8 % with 40 windows and 10 stations in the tests, and more
  for the Rayleigh processing (factors measured with 60 windows and 10
  stations: 0.81 for the vertical component, 0.76 for ARTBF and 0.61
  for ARDS). Compare
  powers with Bartlett or with many windows.
- The ellipticity of the angular processing is distorted when the
  incoherent noise differs between vertical and horizontal components;
  the correction of Wathelet (2024, Section 3.6) is not implemented.
- The empirical uncertainty of the power is available for Bartlett and
  Capon on single components (not MUSIC, not the Rayleigh processing);
  the coordinate uncertainty for Bartlett on single components only.
  The maxima have no uncertainty of their position yet.
- The planning experiments model one plane wave in white noise with
  independent realizations; correlated windows, coloured noise,
  several sources, sampled coordinate errors and the Rayleigh
  processing are not modelled.
- The FDD gives frequencies limited by the analysis step and the
  window resolution; the enhanced FDD gives frequency and damping of
  well separated, lightly damped modes, but its damping estimates
  scatter by about 10 % with 20 minutes of data and are biased by
  leakage (short windows) or by a truncated bell (high `mac`). Close
  modes and harmonics (Jacobsen, Andersen and Brincker, 2007) are
  not treated.

Planned extensions:

- Uncertainty of the position of the maxima (wavenumber, azimuth,
  velocity, ellipticity), to be used also as kernel widths of the
  density; uncertainty of RTBF.
- Correction of the vertical-horizontal noise difference of the
  angular processing (Wathelet, 2024, Section 3.6).
- Bias correction of Capon, $`(L - n + 1) / L`$, as an option.
- Detection of harmonics in the FDD (Jacobsen, Andersen and
  Brincker, 2007); uncertainty of the enhanced FDD over blocks.
- Multitaper spectral estimates with a rigorous number of degrees of
  freedom; the `CrossSpectra` of the FDD used also by the beamforming
  and by spatial autocorrelation (SPAC).
- Mode tracking across frequencies (manual, automatic, hybrid) from
  the density peaks.
- Planning with correlated windows, coloured noise, several sources,
  sampled coordinate errors and RTBF; recommendations on aperture and
  spacing; a constrained geometry optimizer.
- Interactive graphics.

## 20. References

Allemang, R. J., and Brown, D. L. (1982). A correlation coefficient for
modal vector analysis. Proceedings of the 1st International Modal
Analysis Conference, 110-116.

Brincker, R., Ventura, C. E., and Andersen, P. (2001). Damping
estimation by frequency domain decomposition. Proceedings of the 19th
International Modal Analysis Conference, 698-703.

Brincker, R., Zhang, L., and Andersen, P. (2001). Modal identification
of output-only systems using frequency domain decomposition. Smart
Materials and Structures, 10(3), 441-445.

Capon, J. (1969). High-resolution frequency-wavenumber spectrum
analysis. Proceedings of the IEEE, 57(8), 1408-1418.

Capon, J., and Goodman, N. R. (1970). Probability distributions for
estimators of the frequency-wavenumber spectrum. Proceedings of the
IEEE, 58(10), 1785-1786.

Ermert, L., Poggi, V., Burjánek, J., and Fäh, D. (2014). Fundamental
and higher two-dimensional resonance modes of an Alpine valley.
Geophysical Journal International, 198(2), 795-811.

Jacobsen, N.-J., Andersen, P., and Brincker, R. (2007). Eliminating
the influence of harmonic components in operational modal analysis.
Proceedings of the 25th International Modal Analysis Conference
(IMAC).

Percival, D. B., and Walden, A. T. (1993). Spectral Analysis for
Physical Applications. Cambridge University Press.

Schmidt, R. O. (1986). Multiple emitter location and signal parameter
estimation. IEEE Transactions on Antennas and Propagation, 34(3),
276-280.

Wathelet, M. (2024). Incoherent noise-induced distortions of Rayleigh
wave ellipticity measurements obtained with three-component
beamforming. Geophysical Journal International, 236(3), 1804-1827.

Wathelet, M., Guillier, B., Roux, P., Cornou, C., and Ohrnberger, M.
(2018). Rayleigh wave three-component beamforming: signed ellipticity
assessment from high-resolution frequency-wavenumber processing of
ambient vibration arrays. Geophysical Journal International, 215(1),
507-523.

Welch, P. D. (1967). The use of fast Fourier transform for the
estimation of power spectra: a method based on time averaging over
short, modified periodograms. IEEE Transactions on Audio and
Electroacoustics, 15(2), 70-73.
