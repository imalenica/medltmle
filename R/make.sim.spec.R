#' make.sim.spec
#'
#' This function creates functions for iterative expectations.
#'
#' @param timepoint Specifies the final time point.
#' @param YisL Logical indicating whether final outcome is a function of the time-varying covariate.
#'
#' @return Returns a list with models for each time point and covariate.
#'
#' @export make.sim.spec

make.sim.spec <- function(timepoint, YisL=FALSE){

  if(YisL){

    spec <- list(
      qL.c=c(LA_1='LA_1~W1+W2+A_1',LZ_1='LZ_1~W2+A_1+Z_1',Y_1='Y_1~W1+W2+A_1+Z_1+LZ_1+LA_1+A_1:Z_1'),
      qz.c= c(Z_1='Z_1~W2+A_1+LA_1'),
      g.c=c(C_1='C_1~W2+W1',A_1='A_1~W1+W2'),

      QZ.c=c(Z_1='Q.kplus1~W1+W2+A_1+LA_1'),
      QL.c=c(LA_1='Q.kplus1~W1+W2+A_1',LZ_1='Q.kplus1~W1+W2+A_1+LA_1+Z_1+A_1:Z_1'),

      qL.m = c(LA_1='LA_1~A_1',LZ_1='LZ_1~A_1+Z_1',Y_1='Y_1~A_1+Z_1'),
      qz.m = c(Z_1='Z_1~A_1'),
      g.m=c(C_1='C_1~1',A_1='A_1~1'),
      QZ.m =c(Z_1='Q.kplus1~A_1'),
      QL.m =c(LA_1='Q.kplus1~A_1',LZ_1='Q.kplus1~A_1+Z_1'))

    if(timepoint >1){
      for(tt in 2:timepoint){
        spec$qL.c <- c(spec$qL.c, paste0('LA_',tt,'~A_',tt,'+LA_',tt-1,'+LZ_',tt-1),
                       paste0('LZ_',tt,'~W2+A_',tt,'+Z_',tt,'+LZ_',tt-1), paste0('Y_',tt,'~W1+W2+A_',tt,'+Z_',tt,'+LA_',tt,'+LZ_',tt,'+A_',tt,':Z_',tt,'+LA_',tt-1))
        names(spec$qL.c)[length(spec$qL.c)] <- paste0('Y_',tt)
        names(spec$qL.c)[length(spec$qL.c)-1] <- paste0('LZ_',tt)
        names(spec$qL.c)[length(spec$qL.c)-2] <- paste0('LA_',tt)

        spec$qz.c = c(spec$qz.c, paste0('Z_',tt,'~W2+A_',tt,'+LA_',tt))
        names(spec$qz.c)[length(spec$qz.c)] <- paste0('Z_',tt)

        spec$g.c <- c(spec$g.c,paste0('C_',tt,'~W2+A_',tt-1,'+LZ_',tt-1),paste0('A_',tt,'~W2+A_',tt-1,'+LZ_',tt-1))
        names(spec$g.c)[length(spec$g.c)] <- paste0('A_',tt)
        names(spec$g.c)[length(spec$g.c)-1] <- paste0('C_',tt)

        spec$QZ.c <- c(spec$QZ.c, paste0('Q.kplus1~W1+W2+',paste(paste(rep(c('A','Z','LA','LZ'),length(tt-1)),rep(c(1:(tt-1)),each=4),sep='_'),collapse='+'),'+',
                                         paste(paste(paste0('A_',1:(tt-1)),paste0('Z_',1:(tt-1)),sep=':'),collapse='+'),'+A_',tt,'+LA_',tt))
        names(spec$QZ.c)[length(spec$QZ.c)] <- paste0('Z_',tt)

        spec$QL.c <- c(spec$QL.c, paste0('Q.kplus1~W1+W2+',paste(paste(rep(c('A','Z','LA','LZ'),length(tt-1)),rep(c(1:(tt-1)),each=4),sep='_'),collapse='+'),'+',
                                         paste(paste(paste0('A_',1:(tt-1)),paste0('Z_',1:(tt-1)),sep=':'),collapse='+'),'+A_',tt),
                       paste0('Q.kplus1~W1+W2+',paste(paste(rep(c('A','Z','LA','LZ'),length(tt-1)),rep(c(1:(tt-1)),each=4),sep='_'),collapse='+'),'+',
                              paste(paste(paste0('A_',1:(tt)),paste0('Z_',1:(tt)),sep=':'),collapse='+'),'+A_',tt,'+LA_',tt,'+Z_',tt))
        names(spec$QL.c)[length(spec$QL.c)] <- paste0('LZ_',tt)
        names(spec$QL.c)[length(spec$QL.c)-1] <- paste0('LA_',tt)

        spec$qL.m <- c(spec$qL.m, paste0('LA_',tt,'~A_',tt), paste0('LZ_',tt,'~A_',tt,'+Z_',tt), paste0('Y_',tt,'~A_',tt,'+Z_',tt))
        names(spec$qL.m)[length(spec$qL.m)] <- paste0('Y_',tt)
        names(spec$qL.m)[length(spec$qL.m)-1] <- paste0('LZ_',tt)
        names(spec$qL.m)[length(spec$qL.m)-1] <- paste0('LA_',tt)

        spec$qz.m <- c(spec$qz.m, paste0('Z_',tt,'~A_',tt))
        names(spec$qz.m)[length(spec$qz.m)] <- paste0('Z_',tt)

        spec$g.m <- c(spec$g.m,paste0('C_',tt,'~1'),paste0('A_',tt,'~1'))
        names(spec$g.m)[length(spec$g.m)] <- paste0('A_',tt)
        names(spec$g.m)[length(spec$g.m)-1] <- paste0('C_',tt)

        spec$QZ.m <- c(spec$QZ.m, paste0('Q.kplus1~A_',tt))
        names(spec$QZ.m)[length(spec$QZ.m)] <- paste0('Z_',tt)

        spec$QL.m <- c(spec$QL.m, paste0('Q.kplus1~A_',tt), paste0('Q.kplus1~A_',tt,'+Z_',tt))
        names(spec$QL.m)[length(spec$QL.m)] <- paste0('LZ_',tt)
        names(spec$QL.m)[length(spec$QL.m)-1] <- paste0('LA_',tt)
      }
    }

  }else{

    if(timepoint==1){

      spec <- list(
        qL.c=c(LA_1='LA_1~W1+W2+A_1',LZ_1='LZ_1~W2+A_1+Z_1',Y_1='Y_1~W1+W2+A_1+Z_1+LZ_1+LA_1+A_1:Z_1'),
        qz.c= c(Z_1='Z_1~W2+A_1+LA_1'),
        g.c=c(C_1='C_1~W2+W1',A_1='A_1~W1+W2'),

        QZ.c=c(Z_1='Q.kplus1~W1+W2+A_1+LA_1'),
        QL.c=c(LA_1='Q.kplus1~W1+W2+A_1',LZ_1='Q.kplus1~W1+W2+A_1+LA_1+Z_1+A_1:Z_1',Y_1='Q.kplus1~W1+W2+A_1+LA_1+Z_1+A_1:Z_1+LZ_1' ),

        qL.m = c(LA_1='LA_1~A_1',LZ_1='LZ_1~A_1+Z_1',Y_1='Y_1~A_1+Z_1'),
        qz.m = c(Z_1='Z_1~A_1'),
        g.m=c(C_1='C_1~1',A_1='A_1~1'),
        QZ.m =c(Z_1='Q.kplus1~A_1'),
        QL.m =c(LA_1='Q.kplus1~A_1',LZ_1='Q.kplus1~A_1+Z_1'))

    }else{

      spec <- list(
        qL.c=c(LA_1='LA_1~W1+W2+A_1',LZ_1='LZ_1~W2+A_1+Z_1',Y_1='Y_1~W1+W2+A_1+Z_1+LZ_1+LA_1+A_1:Z_1'),
        qz.c= c(Z_1='Z_1~W2+A_1+LA_1'),
        g.c=c(C_1='C_1~W2+W1',A_1='A_1~W1+W2'),

        QZ.c=c(Z_1='Q.kplus1~W1+W2+A_1+LA_1'),
        QL.c=c(LA_1='Q.kplus1~W1+W2+A_1',LZ_1='Q.kplus1~W1+W2+A_1+LA_1+Z_1+A_1:Z_1'),

        qL.m = c(LA_1='LA_1~A_1',LZ_1='LZ_1~A_1+Z_1',Y_1='Y_1~A_1+Z_1'),
        qz.m = c(Z_1='Z_1~A_1'),
        g.m=c(C_1='C_1~1',A_1='A_1~1'),
        QZ.m =c(Z_1='Q.kplus1~A_1'),
        QL.m =c(LA_1='Q.kplus1~A_1',LZ_1='Q.kplus1~A_1+Z_1'))

        for(tt in 2:timepoint){
          spec$qL.c <- c(spec$qL.c, paste0('LA_',tt,'~A_',tt,'+LA_',tt-1,'+LZ_',tt-1),
                         paste0('LZ_',tt,'~W2+A_',tt,'+Z_',tt,'+LZ_',tt-1), paste0('Y_',tt,'~W1+W2+A_',tt,'+Z_',tt,'+LA_',tt,'+LZ_',tt,'+A_',tt,':Z_',tt,'+LA_',tt-1))
          names(spec$qL.c)[length(spec$qL.c)] <- paste0('Y_',tt)
          names(spec$qL.c)[length(spec$qL.c)-1] <- paste0('LZ_',tt)
          names(spec$qL.c)[length(spec$qL.c)-2] <- paste0('LA_',tt)

          spec$qz.c = c(spec$qz.c, paste0('Z_',tt,'~W2+A_',tt,'+LA_',tt))
          names(spec$qz.c)[length(spec$qz.c)] <- paste0('Z_',tt)

          spec$g.c <- c(spec$g.c,paste0('C_',tt,'~W2+A_',tt-1,'+LZ_',tt-1),paste0('A_',tt,'~W2+A_',tt-1,'+LZ_',tt-1))
          names(spec$g.c)[length(spec$g.c)] <- paste0('A_',tt)
          names(spec$g.c)[length(spec$g.c)-1] <- paste0('C_',tt)

          spec$QZ.c <- c(spec$QZ.c, paste0('Q.kplus1~W1+W2+',paste(paste(rep(c('A','Z','LA','LZ'),length(tt-1)),rep(c(1:(tt-1)),each=4),sep='_'),collapse='+'),'+',
                                           paste(paste(paste0('A_',1:(tt-1)),paste0('Z_',1:(tt-1)),sep=':'),collapse='+'),'+A_',tt,'+LA_',tt))
          names(spec$QZ.c)[length(spec$QZ.c)] <- paste0('Z_',tt)

          spec$QL.c <- c(spec$QL.c, paste0('Q.kplus1~W1+W2+',paste(paste(rep(c('A','Z','LA','LZ'),length(tt-1)),rep(c(1:(tt-1)),each=4),sep='_'),collapse='+'),'+',
                                           paste(paste(paste0('A_',1:(tt-1)),paste0('Z_',1:(tt-1)),sep=':'),collapse='+'),'+A_',tt),
                         paste0('Q.kplus1~W1+W2+',paste(paste(rep(c('A','Z','LA','LZ'),length(tt-1)),rep(c(1:(tt-1)),each=4),sep='_'),collapse='+'),'+',
                                paste(paste(paste0('A_',1:(tt)),paste0('Z_',1:(tt)),sep=':'),collapse='+'),'+A_',tt,'+LA_',tt,'+Z_',tt))

          names(spec$QL.c)[length(spec$QL.c)] <- paste0('LZ_',tt)
          names(spec$QL.c)[length(spec$QL.c)-1] <- paste0('LA_',tt)

          spec$qL.m <- c(spec$qL.m, paste0('LA_',tt,'~A_',tt), paste0('LZ_',tt,'~A_',tt,'+Z_',tt), paste0('Y_',tt,'~A_',tt,'+Z_',tt))
          names(spec$qL.m)[length(spec$qL.m)] <- paste0('Y_',tt)
          names(spec$qL.m)[length(spec$qL.m)-1] <- paste0('LZ_',tt)
          names(spec$qL.m)[length(spec$qL.m)-1] <- paste0('LA_',tt)

          spec$qz.m <- c(spec$qz.m, paste0('Z_',tt,'~A_',tt))
          names(spec$qz.m)[length(spec$qz.m)] <- paste0('Z_',tt)

          spec$g.m <- c(spec$g.m,paste0('C_',tt,'~1'),paste0('A_',tt,'~1'))
          names(spec$g.m)[length(spec$g.m)] <- paste0('A_',tt)
          names(spec$g.m)[length(spec$g.m)-1] <- paste0('C_',tt)

          spec$QZ.m <- c(spec$QZ.m, paste0('Q.kplus1~A_',tt))
          names(spec$QZ.m)[length(spec$QZ.m)] <- paste0('Z_',tt)

          spec$QL.m <- c(spec$QL.m, paste0('Q.kplus1~A_',tt), paste0('Q.kplus1~A_',tt,'+Z_',tt))
          names(spec$QL.m)[length(spec$QL.m)] <- paste0('LZ_',tt)
          names(spec$QL.m)[length(spec$QL.m)-1] <- paste0('LA_',tt)
        }

      #Add Y final:
      tt<-timepoint
      spec$QL.c <- c(spec$QL.c, paste0('Q.kplus1~W1+W2+',paste(paste(rep(c('A','Z','LA','LZ'),length(tt-1)),rep(c(1:(tt-1)),each=4),sep='_'),collapse='+'),'+',
                            paste(paste(paste0('A_',1:(tt)),paste0('Z_',1:(tt)),sep=':'),collapse='+'),'+A_',tt,'+LA_',tt,'+Z_',tt,'+LZ_',tt))

      names(spec$QL.c)[length(spec$QL.c)] <- paste0('Y_',tt)


    }

  }

  spec
}

