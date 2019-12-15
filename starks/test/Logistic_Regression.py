import starks
from starks.modp import IntegersModP
from starks.polynomial import polynomials_over
from starks.finitefield import FiniteField
from starks.floatingpoint import FloatingPoint

N = 2

getBin = lambda x: x > 0 and str(bin(x))[2:] or "-" + str(bin(x))[3:]

def floatToBinary64(value):
    val = struct.unpack('Q', struct.pack('d', value))[0]
    return getBin(val).zfill(63)

#reading data from the input file
def Read_from_File (File_Name):
  f = open(File_Name, "r").read().splitlines();

  output_list = []
  temp_list = []

  for line in f:
    num = float(line)
    s = 0
    if num < 0:
      s = 1
      num = num * -1
    num_binary = floatToBinary64(num)

    z = 1
    for j in range(62, 11, -1):
      if num_binary[j] == '1':
        z = 0
      temp_list.append(int(num_binary[j]))
    output_list.append(temp_list)
    temp_list.clear()

    for j in range(10, 0, -1):
      if num_binary[j] == '1':
        z = 0
      temp_list.append(int(num_binary[j]))
    output_list.append(temp_list)
    temp_list.clear()

    temp_list.append(z)
    output_list.append(temp_list)
    temp_list.clear()

    temp_list.append(s)
    output_list.append(temp_list)
    temp_list.clear()

  return output_list
  
# implementation of Logistic Regression in binary finite field
def Logistic_Regression():
  p = 2
  m = 64
  Zp = IntegersModP(p)
  polysOver = polynomials_over(Zp)
  coefficients = [Zp(0)] * 65
  coefficients[0] = Zp(1)
  coefficients[1] = Zp(1)
  coefficients[3] = Zp(1)
  coefficients[4] = Zp(1)
  coefficients[64] = Zp(1)
  poly = polysOver(coefficients)
  field = FiniteField(p, m, polynomialModulus=poly)
  float_number = FloatingPoint(field)

  l = Read_from_File("Input_LR.txt")
  
  w = []
  x = []

  for i in range(N):
    w.append(float_number(field(polysOver(l[4*i])), field(polysOver(l[4*i+1])), l[4*i+2][0], l[4*i+3][0]))

  for i in range(N):
    x.append(float_number(field(polysOver(l[4*N+4*i])), field(polysOver(l[4*N+4*i+1])), l[4*N+4*i+2][0], l[4*N+4*i+3][0]))

  b = float_number(field(polysOver(l[8*N])), field(polysOver(l[8*N+1])), l[8*N+2][0], l[8*N+3][0])
  e = float_number(field(polysOver(l[8*N+4])), field(polysOver(l[8*N+5])), l[8*N+6][0], l[8*N+7][0])


  w_x_mul_add = float_number(field(polysOver([0])), field(polysOver([0])), 0, 0)
  for i in range(N):
    w_x_mul_add = w_x_mul_add + w[i] * x[i]

  w_x_mul_add = w_x_mul_add + b

  #Exp = e ** w_x_mul_add

  #output = Exp/ (Exp + float_number(field(polysOver([1])), field(polysOver([0])), 0, 0))

  #return output
  return w_x_mul_add
