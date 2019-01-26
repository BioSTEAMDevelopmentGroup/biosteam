@echo off 

cd sphinx

start make clean
timeout 5
start make html

Exit