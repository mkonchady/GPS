#*-------------------------------------------------------------------------------
#*- Script to calculate cross correlation of PRN codes for NavIC in the L1 band
#*-------------------------------------------------------------------------------  
import prn_code_navic_l1 as navic
import numpy as np

def octal_to_binary(octal_list):
    binary_values = [bin(int(d, 8))[2:].zfill(3) for d in octal_list]   # convert octal to binary
    result = "".join(binary_values)                                     # make a single string
    return [int(x) for x in list(result)]

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
    
    max_corr = np.max(np.abs(correlation))
    normalized_corr = max_corr / len(p1)
    return normalized_corr * 100.0  # return as percentage


# ----------------------------------------
# 
# ----------------------------------------
if __name__ == "__main__":
    NUM_SATELLITES = 64
    for i in range (0, NUM_SATELLITES):
        tot_corr = 0
        prn1_binary = octal_to_binary(navic.gen_navic_sequence(i+1))
        prn1 = [-1 if x == 0 else x for x in prn1_binary] # convert 0 to -1 for all elements
        zero_pct = prn1.count(-1) / navic.PRN_CODE_LEN 
        for j in range (0, NUM_SATELLITES):
            if (i == j):
                continue   
            prn2_binary= octal_to_binary(navic.gen_navic_sequence(j+1))
            prn2 = [-1 if x == 0 else x for x in prn2_binary] # convert 0 to -1 for all elements
            corr = compare_prn_codes(prn1, prn2)
            tot_corr = tot_corr + corr
        avg_corr = tot_corr / (NUM_SATELLITES - 1)
        print (i, avg_corr, zero_pct)

    print ("Done")
