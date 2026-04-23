import numpy as np

#*------------------------------------------------------------
#*- GPS PRN Code generator using Gold Codes for L1 C/A
#*- Navstar GPS Interface Specification: IS-GPS-200
#*-
#*- https://www.navcen.uscg.gov/sites/default/files/pdf/gps/IS-GPS-200N.pdf
#*-  
#*------------------------------------------------------------

CODE_LEN = 1023       # number of bits in the PRN code

#*-- Gold Codes for Satellite Vehicles
SV = {
   1: [2,6],
   2: [3,7],
   3: [4,8],
   4: [5,9],
   5: [1,9],
   6: [2,10],
   7: [1,8],
   8: [2,9],
   9: [3,10],
  10: [2,3],
  11: [3,4],
  12: [5,6],
  13: [6,7],
  14: [7,8],
  15: [8,9],
  16: [9,10],
  17: [1,4],
  18: [2,5],
  19: [3,6],
  20: [4,7],
  21: [5,8],
  22: [6,9],
  23: [1,3],
  24: [4,6],
  25: [5,7],
  26: [6,8],
  27: [7,9],
  28: [8,10],
  29: [1,6],
  30: [2,7],
  31: [3,8],
  32: [4,9],
}

#*--- add 1 to the list of first octals to get the first 10 chips
FIRST_OCTALS = ['440', '620', '710', '744', '133', '455', '131', '454', '626', '504', '642', '750', '764', 
                '772', '775', '776', '156', '467', '633', '715', '746', '763', '063', '706', '743', '761', 
                '770', '774', '127', '453', '625', '712']

#*-- convert a 3 digit octal to a binary
BINARY_DIGITS = 9
def octal_to_binary(octal_str):
    decimal_num = int(octal_str, 8)         # convert octal base 8 to decimal
    binary_num = bin(decimal_num)[2:]       # convert decimal to binary
    return "0" * (BINARY_DIGITS - len(binary_num)) + binary_num

#*---------------------------------------------------------------------------------------
#*- Parameters:
#*-    register: 10 element array
#*-    feedback: List of positions that are added modulo 2 and pre-pended to position 0
#*-    result:   List of positions that are added module 2 and returned to caller
#*-  Returns the bit that will be suffixed to the PRN code
#*----------------------------------------------------------------------------------------
def shift(register, feedback_pos, result_pos):    
    # calculate the value of the bit that will be added to the PRN code
    res_arr = [register[i-1] for i in result_pos]
    if len(res_arr) > 1:
        result = sum(res_arr)  % 2
    else:
        result = res_arr[0]
        
    # calculate the feedback bit that will be inserted in position 0
    feedback = sum([register[i-1] for i in feedback_pos]) % 2
    
    # shift 1 position to the right
    for i in reversed(range(len(register[1:]))):
        register[i+1] = register[i]
    register[0] = feedback
    
    return result

#*---------------------------------------------------------------------------------------
#*- Parameters:
#*-    sv: Satellite Vehicle number
#*- Returns the PRN Code binary string where 0s are converted to -1s
#*----------------------------------------------------------------------------------------
def PRN(sv):
    
    global FIRST_OCTALS

    # initialize the two registers
    REGISTER_LEN = 10
    G1 = [1 for i in range(REGISTER_LEN)]
    G2 = [1 for i in range(REGISTER_LEN)]

    prn_code = [] # list of PRN code bits
    
    # Generate a bit that will be added to ca
    for i in range(CODE_LEN):
        g1 = shift(G1, [3,10], [10])            # Generate the first part from the shift
        g2 = shift(G2, [2,3,6,8,9,10], SV[sv])  # Generate the second part from the shift for the particular satellite
        prn_code.append((g1 + g2) % 2)
    
    # verify that the first 10 chips of the PRN code match the first 10 chips from the spec
    first_10_chips = '1' + octal_to_binary(FIRST_OCTALS[sv - 1])
    first_10_prncode = "".join([str(x) for x in prn_code[:10]])
    assert first_10_chips == first_10_prncode, "First 10 chips of " + str(sv) + " did not match"

    return  [-1 if x == 0 else x for x in prn_code] # convert 0 to -1 for all elements

#*------------------------------------------------------------------
#*-- Compares two GPS PRN codes using normalized cross-correlation.
#*-- Assumes PRN codes are mapped as 0 -> -1 and 1 -> 1.
#*------------------------------------------------------------------
def compare_prn_codes(prn1, prn2):
    
    p1 = np.array(prn1)
    p2 = np.array(prn2)

    #*-- the product of ffts results in the frequency domain of the cross correlation
    correlation = np.fft.ifft(np.fft.fft(p1) * np.conj(np.fft.fft(p2)))
    correlation = np.real(correlation)
    
    # 3. Find max correlation
    max_corr = np.max(np.abs(correlation))
    
    # 4. Normalize (if codes are 1023 bits, max is 1023)
    normalized_corr = max_corr / len(p1)
    
    return normalized_corr * 100.0  # return as percentage

# Main 

tot_corr = 0
NUM_SATELLITES = 32
for i in range (0, NUM_SATELLITES):
    tot_corr = 0
    for j in range (0, NUM_SATELLITES):
        if (i == j):
            continue
        prn1 = PRN(i+1)
        prn2 = PRN(j+1)

        corr = compare_prn_codes(prn1, prn2)
        tot_corr = tot_corr + corr
    print (i, tot_corr / (NUM_SATELLITES - 1))
