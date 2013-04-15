#  Converted from MASM to GAS
#
#	.	686
#	.	model	flat
#	.	code
.text
      # i=merge(*sx,lx,*sy,ly,*sz, lcol, *ncomp), with updates to ncomp
      #
      # merge annotated string pointers at sx, sy into sz.
      # sz and sy may overlap at the top.
      # counts of items in sx,sy are in lx,ly
      #
      # after setup,
      # eax points to sx; [esp]   is end of sx
      # ecx points to sy; [esp+4] is end of sy
      # edi points to sz

	.globl	_merge            # _cdecl entry point

_merge:	#proc
	mov	4(%esp),%eax         # x  for _cdecl
	mov	8(%esp),%edx         # lx for _cdecl
	mov	12(%esp),%ecx        # y  for _cdecl
	jmp	@merge
#_merge	endp

	.globl	_lcmp             # _cdecl entry point

_lcmp:	#proc
	mov	4(%esp),%eax         # px
	mov	8(%esp),%edx         # py
	jmp	@lcmp
#_lcmp	endp

@merge:	#proc
	push	%edi
	push	%esi
	push	%ebx                 # for optlink
	push	%ebp
                              # set up limits
	sub	$8,%esp              # for 2 pointers
	mov	40(%esp),%ebx        # ly
	shl	$2,%ebx
	add	%ecx,%ebx
	mov	%ebx,(%esp)          # end of sy

	mov	%edx,%ebx            # lx
	shl	$2,%ebx
	add	%eax,%ebx
	mov	%ebx,4(%esp)         # end of sx
	mov	44(%esp),%edi        # sz

@lphd:
	cmp	4(%esp),%eax
	jae	@lab5                # x exhausted
	cmp	(%esp),%ecx
	jae	@lab4                # y exhausted
	mov	(%eax),%ebx          # pointer at x_i
	mov	(%ebx),%ebp          # dereference, getting x.llcp
	mov	(%ecx),%esi          # pointer at y_j
	cmp	(%esi),%ebp          # deference, getting y.llcp
	jg	@lab1                   # x.llcp > y.llcp
	jl	@lab2                   # x.llcp < y.llcp
	push	%eax
	push	%edx
	push	%ecx                 # save around function call
                              # eax, edx parameters to lcmp, optlink linkage
	push	%eax
	mov	68(%esp),%eax        # pointer to ncomp
	incl	(%eax)               # update ncomp
	pop	%eax
	mov	%ecx,%edx
	mov	60(%esp),%ecx        # lcol
	call	@lcmp
	cmp	$0,%eax
	pop	%ecx
	pop	%edx
	pop	%eax
	jg	@lab2                   # not in order, so output y first
@lab1:
	mov	(%eax),%ebx
	add	$4,%eax              # bump up pointer to sx
	jmp	@lab3
@lab2:
	mov	(%ecx),%ebx
	add	$4,%ecx              # bump up pointer to sy
@lab3:
	mov	%ebx,(%edi)          # set z_k to a_i or b_j, as proper
	add	$4,%edi
	jmp	@lphd                # loop
@lab4:
	cmp	4(%esp),%eax         # x exhausted?
	je	@lab5
	mov	(%eax),%ebx
	mov	%ebx,(%edi)          # z <- x
	add	$4,%eax
	add	$4,%edi
	jmp	@lab4
@lab5:
	add	$8,%esp
	pop	%ebp
	pop	%ebx
	pop	%esi
	pop	%edi
	ret
#	@merge	endp
#
      # optlink entry point
      # i=lcmp(*sx,*sy,lcol). sx and sy are {int llcp; char *str;}
      #
      # compare strings starting at (sx.str)+sx.llcp, (sy.str)+sx.llcp
      # stop at null or first non-matching character or at lcol.
      # if first non-matching character of sx is greater, return 1, update sx.llcp
      # otherwise, return -1, update sy.llcp
      # should not be called with sx.llcp different from sy.llcp
      #

@lcmp:	#proc
	push	%edi
	push	%esi
	push	%ebx                 # for optlink
	push	%ebp                 #
	mov	(%eax),%eax          # dereference
	mov	(%edx),%edx
	mov	(%eax),%esi          # sx.llcp
	mov	(%edx),%edi          # sy.llcp
	mov	%edx,%ebx            # save offset of second par
	cmp	%edi,%esi
	jz	@blb1
	xor	%eax,%eax
	jmp	@blbl10
@blb1:
	mov	4(%eax),%ebp         # base of sx.str
	mov	4(%edx),%edi         # sy.str
	add	%esi,%edi            # first character of stry to compare
	lea	(%ebp,%esi),%esi     # sx.str+sx.llcp
	lea	(%ecx,%ebp),%ecx     # beyond lcol
	cld
@next:
	cmp	%ecx,%esi
	jne	@nx                  # at lcol+1 ?
	dec	%esi
	jmp	@done
@nx:
	movb	(%esi),%dl
	orb	%dl,%dl
	jz	@done                   # end of sx.str?
	cmpsb
	je	@next
	dec	%esi                 # to compensate for autoincrement
@done:
	sub	%ebp,%esi            # new llcp of first string
	cmpb	-1(%edi),%dl
	jbe	@less                # first string is less
	mov	%esi,(%eax)          # eax has not been changed
	mov	$1,%eax
	jmp	@blbl10
@less:
	mov	%esi,(%ebx)          # offset of second parameter
	mov	$-1,%eax
@blbl10:
	pop	%ebp
	pop	%ebx
	pop	%esi
	pop	%edi
	ret
#@lcmp	endp
#	end
