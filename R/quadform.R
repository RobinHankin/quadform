`ht` <- function(x){ t(Conj(x)) }

`cprod` <- function(x,y=NULL){
    if(is.null(y)){
        return(crossprod(Conj(x),x))
    } else {
        return(crossprod(Conj(x),y))
    }
}

`tcprod` <- function(x,y=NULL){
    if(is.null(y)){
        return(tcrossprod(x,Conj(x)))
    } else {
        return(tcrossprod(x,Conj(y)))
    }
}

`quad.form.chol` <- function (chol, x){
        jj <- cprod(chol, x)
        drop(cprod(jj, jj))
}

`quad.form` <- function (M, x){ drop(crossprod(crossprod(M,Conj(x)),x)) }

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

`quad.trace` <- function(M,x){ sum(crossprod(M, Conj(x)) * x) }

`quad.ttrace` <- function(M,x){ sum(tcrossprod(Conj(x), M) * x) }

cp <- cprod
tcp <- tcprod
qf <- quad.form
qfi <- quad.form.inv
q3 <- quad.3form
q3i <- quad.3form.inv

q3t <- quad.3tform
qt <- quad.tform
q3i <- quad.tform.inv
qd <- quad.diag
qtd <- quad.tdiag
q3d <- quad.3diag
q3td <- quad.3tdiag

qtr <- quad.trace
qttr <- quad.ttrace
