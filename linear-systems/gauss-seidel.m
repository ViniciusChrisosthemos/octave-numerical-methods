## Copyright (C) 2018, Anderson R. P. Domingues
##
## This file (seidel.m) is free software; you can redistribute 
## it and/or modify it under the terms of the GNU General Public
## License as published by the Free Software Foundation;
## either version 3 of the License, or (at your option) any
## later version.
##
## This file (seidel.m) is distributed in the hope that it will 
## be useful, but WITHOUT ANY WARRANTY; without even the implied
## warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
## PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public
## License along with Octave; see the file COPYING.  If not,
## see <http://www.gnu.org/licenses/>.

## usage: X = seidel (A, B, it, ex, init)
##
## Solve a linear system using Seidel method. Parameter <A> is 
## coeficient matrix, <B> is the result vector, <it> is the 
## maximum number of iterations (used as stop criteria), <ex> 
## is the number of exact digits (stop criteria), and <init>
## is the partial solution vector (used as partial results for
## the first iteration). One may ommit the <init> param, which
## the default is zeroes(column(a), rows(a)) -- see columns, 
## rows, and zeroes functions.
##
## Example:
##
## A = [
##   10   -1    2    0
##   -1   11   -1    3
##    2   -1   10   -1
##    0    2   -1    8]
##
## B = [ 
##    6
##   25
##  -11
##   15]
##
## X = seidel(A, B, 10, 6)
##
## >> ans = [
##   0.98827
##   1.92704
##  -0.97785
##   1.27101 
## ]

## Author: Anderson R. P. Domingues
## Keywords: linear systems seidel 
## Maintainer: Anderson R. P. Domingues
## Email: anderson.domingues@acad.pucrs.br


##Função que resolve um sistema de equações com o método de seidel
#Gera uma matrix auxiliar, com uma coluna de 1´s a esquerda, zera a diagonal e troca o sinal dos valores
#Gera um vetor [1,a0,b0,...,z0] para multiplicar as linhha da matriz auxiliar
#Repete até i > MAXIT ou ErroTotal > RTOL
function retval = seidel(A, B, MAXIT = 20, RTOL = 1e-6 ,X0 = zeros(1,columns(A)))  
  #Inicia variaveis
  length = columns(A);                                      
  diagonal_a = transpose(diag(A));                          
  values = ones(1,length+1) ;                               
  matrix = ones(length,length+1);           

  #Formata matriz auxiliar  
  for i=1 : length                              
    values(i+1) = X0(i);                                    
                                                            
    matrix(i,1) = B(i,1);                                   
    matrix(:,i+1) = [A(:,i)*-1];                            
    matrix(i,i+1) = 0;                                      
  endfor
  
  #Formata vetor de valores
  values = transpose(values);
  erro = Inf;
  temp = ones(1,length);
  
  it = 1;
  #Laço que calcula an,bn,...,zn
  while (it < MAXIT) && erro > RTOL
    erro = 0;
    
    for value=2 : length+1;
      values(value) = (matrix(value-1,:) * values)/diagonal_a(value-1);
    endfor
    it++;
    
    #Calcula Erro total
    for j=1 : length
      erro += (temp * A(:,j)) * values(j+1);
    endfor
    erro -= temp * B;
    erro = abs(erro);
    
  endwhile
  
  result = zeros(length,1);
  
  #Formata resultado 
  for v=1: length;
    result(v,1) = values(v+1);
  endfor
  
  retval = result
endfunction



if(__SEIDEL_METHOD_TEST == 1)
  %entrada da matriz de coeficientes
  A = [
     10   -1    2    0
     -1   11   -1    3
      2   -1   10   -1
      0    2   -1    8]

  %entrada do vetor solucao    
  B = [ 
      6
     25
    -11
     15]
   
  %funcao para o metodo de Seidel
  X = seidel(A, B);

endif