
 # Mastering-Parallel-Programming-with-R

You can read more at [Mastering Parallel Programming with R]
(https://www.packtpub.com/big-data-and-business-intelligence/mastering-parallel-programming-r?utm_source=Github&utm_medium=Repository&utm_campaign=9781784394004)

# Instructions
#
## Copyright 2016 Simon Chapple
##
##Packt: "Mastering Parallelism with R"
##Chapter 1 - Aristotle's Number Puzzle
##

segueSC replacement for segue
=============================

Unfortunately the original segue package no longer works with the current AWS Java API.
There have also been various changes that have been made to the EMR Hadoop task launcher,
cluster configuration tools and options, AWS security and S3 usage which require wholesale
changes to be made to the implementation of segue.

So dear reader, after a lot of effort, I now present you with my modified version of the
original segue package that does indeed work with today's AWS. This has been most recently
tested in full on 22-May-2016.

segueSC is a drop in replacement for segue.
You can install it directly from this source code tree with:
R> install.packages("<your directory path>/segueSC",repos=NULL,type="source")
R> library("segueSC")

PLEASE remember that you should always check to see what resources you have left running
from your AWS web console. This includes EMR clusters as well as S3 buckets. If you leave
these running then you will be charged for their use. You can always terminate clusters
through the AWS web console manually, and can also delete S3 buckets manually. An EMR
cluster that is terminated will not incur further charges, but any related S3 buckets
that remain will incur charges. So PLEASE DO CHECK whenever you finish your AWS work.

Whilst I have not gone into these aspects in the book itself (since it detracts too much
from the core theme of parallel programming), its certainly instructional to look through
the implementation of the segue API in awsFunctions.R to see how the AWS API is used to
enable remote creation, execution and control of AWS Elastic Map Reduce Hadoop clusters.
 
Some things its useful to know about the segueSC EMR clusters:
AMI version: 3.11.0
Hadoop distribution: Amazon 2.4.0
AWS region used: us-east-1c
Logging is turned on by default.

Simon R Chapple
22-May-2016

If you wish to contact me, you can find me on LinkedIn:
https://uk.linkedin.com/in/simonchapple



Related R Books and Videos

* [R for Data Science] 
* (https://www.packtpub.com/big-data-and-business-intelligence/r-data-science?utm_source=Github&utm_medium=Repository&utm_campaign=9781784390860)
* [Learning R for Data Visualization  [Video]] (https://www.packtpub.com/big-data-and-business-intelligence/learning-r-data-visualization-video?utm_source=Github&utm_medium=Repository&utm_campaign=9781785882890)
* [Machine Learning with R Cookbook] (https://www.packtpub.com/big-data-and-business-intelligence/machine-learning-r-cookbook?utm_source=Github&utm_medium=Repository&utm_campaign=9781783982042)


