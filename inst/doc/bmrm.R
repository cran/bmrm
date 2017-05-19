## ---- echo=TRUE, results='hide'------------------------------------------
    library(bmrm)
    x <- cbind(intercept=100,data.matrix(iris[c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]))
    y <- iris$Species
    w <- nrbm(softMarginVectorLoss(x,iris$Species))
    w <- matrix(w,ncol(x),dimnames=list(colnames(x),levels(iris$Species)))
    predictions <- colnames(w)[max.col(x %*% w)]
    table(target=y,prediction=predictions)

