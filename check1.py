from mpmath import mp, mpf, log, exp, sqrt, pi, cosh, euler, quad, inf

mp.dps = 30                          # working digits; raise to match a weil_step1.py run if you like

lam = mpf(2)
L = 2*log(lam)                       # window length in y = log u

def phi(y):                          # f = V_0 (constant): autocorrelation is a triangle
    return 2*(1 - y/L) if y < L else mpf(0)

def integrand(y):                    # archimedean integrand from section 2 of the note
    return (exp(y/2)*phi(y) - phi(0)) / (exp(y) - exp(-y))

# --- archimedean term: boundary piece + integral over [0,L] + integral over [L,inf) ---
boundary = (euler + log(4*pi))/2 * phi(0)
part1 = quad(integrand, [0, L])      # quad handles the removable 0/0 at y=0 (tanh-sinh never evaluates the endpoint)
part2 = quad(integrand, [L, inf])    # phi = 0 here; mpmath integrates to infinity directly
arch = boundary + part1 + part2

# --- pole term (closed form) and prime term (prime powers <= lambda^2 = 4) ---
pole = 16/L*(cosh(L/2) - 1)
primes = (log(2)/sqrt(2))*phi(log(2)) + (log(3)/sqrt(3))*phi(log(3)) + (log(2)/2)*phi(2*log(2))

print("archimedean", arch)           # expect 2.0623050824749...
print("  tail check", part2, "vs closed form", -log(mpf(5)/3))   # integral over [L,inf) = -ln coth(L/2) = -ln(5/3) at lambda=2
print("pole       ", pole)           # expect 2.8853900817779...
print("primes     ", primes)         # expect 0.7533807586786...
print("a0 = pole - arch - primes =", pole - arch - primes)       # expect 0.069704243... -- compare to 'coefficients n= 0  a=' in a lambda=2 Fourier log