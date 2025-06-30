import sys
import os
import random

import pydivsufsort
import numpy as np

def sim_seq(length):
	alphabet = {0:"A", 1:"C", 2:"G", 3:"T"}
	seq = []
	for r in range(0, length):
		randchar = random.randint(0, 3)
		seq.append(alphabet[randchar])
	seq = ''.join(seq)	
	
	return(seq)

def tsd_err(seq, errprob = 1, max_err = 2):
	alphabet = {0:"A", 1:"C", 2:"G", 3:"T"}
	seq = list(seq)
	newseq = []
	for c in seq:
		pull = random.randint(0, 10000)
		if pull < errprob:
			newchar = c
			while newchar == c:
				new_pull = random.randint(0, 3)
				newchar = alphabet[new_pull]
			
			c = newchar
		newseq.append(c)
		
	newseq = ''.join(newseq)
	return newseq


def simulate_tsd(min_tsd_length = 6, max_tsd_length = 250, max_upstream = 15, max_downstream = 15, line_min_length = 1000, line_max_length = 10000, crap_left = 1000, crap_right = 1000):
	tsd_length = random.randint(min_tsd_length, max_tsd_length)
	line_length = random.randint(line_min_length, line_max_length)
	upstream = random.randint(0, max_upstream)
	downstream = random.randint(0, max_downstream)
	crap_l =	random.randint(0, crap_left)
	crap_r = 	random.randint(0, crap_right)
	
	
	#Simulate end finding with some error
	line_end = tsd_length + upstream + line_length + crap_l + random.randint(-3, 3)
	
	line_element = sim_seq(line_length)
	tsd = sim_seq(tsd_length)
	up = sim_seq(upstream)
	down = sim_seq(downstream)
	cl = sim_seq(crap_l)
	cr = sim_seq(crap_r)
	
	left_tsd_err = tsd_err(tsd)
	
	combined_element = ''.join([cl, left_tsd_err, up, line_element, down, tsd, cr])
	
	return left_tsd_err, tsd, line_end, combined_element
	

#Here is where we develop functions for TSD extraction
	
import numpy as np

from collections import defaultdict
#Deepseek's suggested solution and it's pretty good.

#Augmenting this with a method to join adjacent repeats with no more than an N-character mismatch, 
#one insertion followed by continued matching would probably be good.
#Should pass this function a string pasted together without the middle line element involved, since we find the ends of it either way.

def compute_suffix_array(s):
    n = len(s)
    SA = sorted(range(n), key=lambda i: s[i:])
    return SA

def compute_lcp_array(s, SA):
    n = len(s)
    LCP = [0] * n
    for i in range(1, n):
        a, b = SA[i-1], SA[i]
        h = 0
        while a + h < n and b + h < n and s[a + h] == s[b + h]:
            h += 1
        LCP[i] = h
    return LCP

def find_approximate_repeats(s, min_L=5, mismatch_max = 2):
    n = len(s)
    if n == 0:
        return [], []
    SA = compute_suffix_array(s)
    LCP = compute_lcp_array(s, SA)
    resultsA = []
    resultsB = []
    
    for i in range(1, n):
        p1, p2 = SA[i-1], SA[i]
        h = LCP[i]
        mismatches = 0
        j = h
        max_j = min(n - p1, n - p2)
		
        while j < max_j and mismatches <= mismatch_max:
            if s[p1 + j] != s[p2 + j]:
                mismatches += 1
                if mismatches > mismatch_max:
                    break
            j += 1
        L_length = j
        if L_length >= min_L:
            resultsA.append((p1, p2, L_length))
        
        for d in range(1, 11):
            j1 = 0
            while (p1 + h + d + j1 < n and 
                   p2 + h + j1 < n and 
                   s[p1 + h + d + j1] == s[p2 + h + j1]):
                j1 += 1
            if j1 > 0:
                L1 = h + d + j1
                L2 = h + j1
                if L1 >= min_L and L2 >= min_L:
                    resultsB.append((p1, p2, L1, L2))
                    
            j2 = 0
            while (p1 + h + j2 < n and 
                   p2 + h + d + j2 < n and 
                   s[p1 + h + j2] == s[p2 + h + d + j2]):
                j2 += 1
            if j2 > 0:
                L1 = h + j2
                L2 = h + d + j2
                if L1 >= min_L and L2 >= min_L:
                    resultsB.append((p1, p2, L1, L2))
                    
    return resultsA, resultsB

'''
# Example usage
if __name__ == "__main__":
    s = "ABCDEFABCXYZDEF"
    min_L = 5
    resultsA, resultsB = find_approximate_repeats(s, min_L)
    print("Condition A (mismatches):")
    for p1, p2, L in resultsA:
        print(f"Repeated substring of length {L} at positions {p1} and {p2}:")
        print(f"  s1: {s[p1:p1+L]}")
        print(f"  s2: {s[p2:p2+L]}")
    print("\nCondition B (indels):")
    for p1, p2, L1, L2 in resultsB:
        print(f"Repeated substring with lengths {L1} and {L2} at positions {p1} and {p2}:")
        print(f"  s1: {s[p1:p1+L1]}")
        print(f"  s2: {s[p2:p2+L2]}")
'''


class Solution:
	def find_repeated_substrings(self, s, min_length=2):
		n = len(s)
		if n == 0:
			return {}
		
		sa = pydivsufsort.divsufsort(s)
		lcp = pydivsufsort.kasai(s, sa)
		
		max_rep = [0] * n
		max_rep = np.zeros(shape = (n), dtype = np.int32)
		
		if n > 1:
			max_rep[sa[0]] = lcp[0]
			max_rep[sa[n-1]] = lcp[n-2]
		
		#Convert loops to vector ops for efficiency
		#For each suffix in the lex. sorted suffix array, the max of adjacent values in the least common prefix array is the largest repeat; 
		#if this value is 0, the string does not repeat anywhere
		lcp_lower = lcp[:-2]
		lcp_upper = lcp[1:-1]
		lcp_maxes = np.maximum(lcp_lower, lcp_upper)
		
		#Longest repeat at each position in the array
		max_rep[sa[1:-1]] = lcp_maxes
	
		
	
		resultsA = []
	
		#This code technically has a limitation in that it only finds adjacent repeats. 
		#There are possible cases where repeats may have multiple mismatches spread out over a significant distance that this will miss.
		#Can we increase 
		
		mismatch_pct = 0.05
		for i in range(1, n):
			p1, p2 = sa[i-1], sa[i]
			h = lcp[i]
			mismatches = 0
			j = h
			max_j = min(n - p1, n - p2)
			
			mismatch_max_pct = int((mismatch_pct * max_j) + 0.5)
			mismatch_max_num = 2
			mismatch_max = min([mismatch_max_pct, mismatch_max_num])
			
			while j < max_j and mismatches <= mismatch_max:
				if s[p1 + j] != s[p2 + j]:
					mismatches += 1
					if mismatches > mismatch_max:
						break
				j += 1
			L_length = j
			if L_length >= min_length:
				resultsA.append((p1, p2, L_length))
		
		resultsB = []
		for i in range(1, n):
			p1, p2 = sa[i-1], sa[i]
			h = lcp[i]
			for d in range(1, 11):
				j1 = 0
				while (p1 + h + d + j1 < n and 
					   p2 + h + j1 < n and 
					   s[p1 + h + d + j1] == s[p2 + h + j1]):
					j1 += 1
				if j1 > 0:
					L1 = h + d + j1
					L2 = h + j1
					if L1 >= min_length and L2 >= min_length:
						resultsB.append((p1, p2, L1, L2))
				
				j2 = 0
				while (p1 + h + j2 < n and 
					   p2 + h + d + j2 < n and 
					   s[p1 + h + j2] == s[p2 + h + d + j2]):
					j2 += 1
				if j2 > 0:
					L1 = h + j2
					L2 = h + d + j2
					if L1 >= min_length and L2 >= min_length:
						resultsB.append((p1, p2, L1, L2))
			
		
		'''
		substr_dict = defaultdict(list)
		for i in range(n):
			if max_rep[i] > 0:
				substr = s[i:i+max_rep[i]]
				substr_dict[substr].append(i)
				
		return dict(substr_dict)
		'''
		
		return resultsA, resultsB
	
mn = Solution()

for i in range(0, 2000):
	true_left_tsd, true_right_tsd, threeprime_line_bound, full_line  = simulate_tsd()
	
	res_a, res_b = mn.find_repeated_substrings(full_line, min_length = 6)
	
	for a in res_a:
		print(a)
	
	longest_substring = 0
	
	biggest_recovery = None
	
	'''
	for d in rd:
		if len(d) > longest_substring:
			longest_substring = len(d)
			biggest_recovery = d
	
	if true_left_tsd != true_right_tsd:
		#print("tsd_mutation", d, len(d), len(true_right_tsd))
		pass
	else:
		#print("tsd_duplicate", d == true_right_tsd, len(d), len(true_right_tsd))
		if d != true_right_tsd:
			#print(d)
			#print(true_right_tsd)
			#print("")
			pass
	'''		
	
	#search_zone = full_line[(threeprime_line_bound-5):]
	
	