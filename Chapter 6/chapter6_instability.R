##
## Copyright 2015 Simon Chapple
##
## Packt: "Mastering Parallelism with R"
## Chapter 6 - Numerical Instability
##

# Create a long vector of fractional numbers
v <- 1:500000
for (i in 1:length(v))
{
  v[i] = 1/i
}

# sum forward
suma <- 0.0
for (i in 1:length(v))
{
  suma = suma + v[i]
}
print(suma,digits=15)

# sum reverse
sumz <- 0.0
for (i in length(v):1)
{
  sumz = sumz + v[i]
}
print(sumz,digits=15)

##
## Packt: "Mastering Parallelism with R"
## Chapter 7 - Numerical Instability
##
## Copyright 2015 Simon Chapple
##
