# the dimensional ODE model written for deSolve

m_cfs <- function(t,y,p)
{
  with( as.list(c( y, p )), {
    
    return(list(c(
      
      dC_dt <- a * C * (P/b - 1) - e * C, 
      
      dI_dt <- f*g*C - w*I - I/(D*s + I)*D + k*(h - f*g*C),
      
      dD_dt <- m*(h*q/P - D),
      
      dP_dt <- r*P*( s*D/I - 1)
      
    )))
    
  })
}

# the non-dimensionalised ODE model written for deSolve
m_cfs_nd <- function(t,y,p)
{
  with( as.list(c( y, p )), {
    
    
    return(list(c(
      dChat_dthat <- Chat * (alpha*Phat- 1) - beta * Chat, 
      
      dIhat_dthat <- delta*Chat - omega*Ihat - gamma*Ihat*Dhat/(Dhat + Ihat) + kappa*(gamma - delta*Chat) , 
      
      dDhat_dthat <- mu*(1/Phat - Dhat),
      
      dPhat_dthat <- rho*Phat*( Dhat/Ihat - 1)
    )))
    
    
  })
}