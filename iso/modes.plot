MACRO pspec 3	{ # Plot spectrum from history points
		# Arg 1 = file to read from
		# Arg 2 = offset from start of file
		# Arg 3 = length of individual spectrum
		define _start ($2)
		define _end   ($2 + $3)
		data $1
		lines $_start $_end
		read { _m 2 _e 3 }
		set _lm = lg(_m)
		set _le = lg(_e)
		connect _lm _le
		foreach v ( _start _end _m _lm _e _le ) { delete $v }
		}

erase
expand 1.5
limits 0 1.21 -8 -1
NOTATION 1 1 1 1
LOCATION 5000 31000 5000 31000
ticksize -1 0 -1 0
box
xlabel \i k
ylabel \i E(k)

do k = 0, 19 {
	define s (3 + 15 * $k)
	pspec iso.his $s 14
}

