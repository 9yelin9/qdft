import re

def re_int(pattern, string):
	return int(re.sub(pattern, '', re.search(f'{pattern}[-]?\d+', string).group()))

def re_float(pattern, string):
	return float(re.sub(pattern, '', re.search(f'{pattern}[-]?\d+[.]\d+', string).group()))
