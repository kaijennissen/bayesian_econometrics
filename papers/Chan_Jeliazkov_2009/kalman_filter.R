library(Matrix)


# Variante 1
# y(t) = F * x(t) + u(t)
# x(t) = G * x_(t-1) + v(t)


plot(AirPassengers)
y <- c(AirPassengers)

# Variante 2

# nu = P^-1 (K* nu_tilde ü G' I)

K <-  t(H)%*% solve(S) %*% H
G <-  bdiag(G_t) # TODO
Omega_11 <-  diag(TT)*1
P <-  K + t(G) %*% kronecker(diag(TT), Omega_11) %*% G
nu_tilde <- 
y <- 
X <- 
beta <- 
nu_hat <- solve(P)*(K*nu_tilde+t(G) %*% kronecker(diag(TT), Omega_11) %*% (y-X%*%beta))
 