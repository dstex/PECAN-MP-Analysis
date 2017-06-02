from numpy import mod,floor

def hhmmss2sec(timehhmmss):
	sec = floor(timehhmmss/10000)*3600 + floor(mod(timehhmmss,10000)/100)*60 + floor(mod(timehhmmss,100))
	
	return sec