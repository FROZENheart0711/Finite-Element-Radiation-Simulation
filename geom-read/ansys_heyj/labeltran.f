	subroutine labeltran(i,j,mpresent)
	if (i.eq.1) then
	mpresent=j-1
	else
	endif
	if (i.eq.2) then
	mpresent=j+1
	else
	endif
	if (i.eq.3) then
	mpresent=6
	else
	endif
	return
	end
