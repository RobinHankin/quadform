`ht` <- function(x){
    if(is.complex(x)){
        return(t(Conj(x)))
    } else {
        return(t(x))
    }
}

`cprod` <- function(x,y=NULL){
    if(is.complex(x) | is.complex(y)){
        if(is.null(y)){
            return(crossprod(Conj(x),x))
        } else {
            return(crossprod(Conj(x),y))
        }
    } else {
        return(crossprod(x,y))
    }
}

`tcprod` <- function(x,y=NULL){
    if(is.complex(x) | is.complex(y)){
        if(is.null(y)){
            return(tcrossprod(x,Conj(x)))
        } else {
            return(tcrossprod(x,Conj(y)))
        }
    } else {
        return(tcrossprod(x,y))
    }
}

`quad.form` <- function (M, x, chol = FALSE){
    if (chol == FALSE) {
        return(drop(crossprod(crossprod(M,Conj(x)),x)))
    }
    else {
        jj <- cprod(M, x)
        return(drop(cprod(jj, jj)))
    }
}

`quad.form.inv` <- function (M, x){ drop(cprod(x, solve(M, x))) }

`quad.3form` <- function(M,left,right){ crossprod(crossprod(M, Conj(left)), right) }

`quad.3form.inv` <- function(M,left,right){ drop(cprod(left, solve(M, right))) }

`quad.3tform` <- function(M,left,right){ tcrossprod(left, tcrossprod(Conj(right), M)) }

`quad.tform` <- function(M,x){ tcrossprod(x, tcrossprod(Conj(x), M)) }

`quad.tform.inv` <- function(M,x){ drop(quad.form.inv(M, ht(x))) }

`quad.diag` <- function(M,x){ colSums(crossprod(M, Conj(x)) * x) }

`quad.tdiag` <- function(M,x){ rowSums(tcrossprod(Conj(x), M) * x) }

`quad.3diag` <- function(M,left,right){ colSums(crossprod(M, Conj(left)) * right) }

`quad.3tdiag` <- function(M,left,right){ colSums(t(left) * tcprod(M, right)) }
