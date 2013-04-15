      .686
      .model flat
      .code
      ; Can be called using OPTLINK or CDECL CONVENTION
      ; i=merge(*sX,lX,*sY,lY,*sZ, lcol, *ncomp), with updates to ncomp
      ;
      ; merge annotated string pointers at sX, sY into sZ.
      ; sZ and sY may overlap at the top.
      ; numbers of items in sX,sY are in lX,lY
      ;
      ; After setup,
      ; eax points to sX; [esp]   is end of sX
      ; ecx points to sY; [esp+4] is end of sY
      ; edi points to sZ

public _merge                              ; _CDECL entry point

_merge proc
      mov     eax,dword ptr[esp+4]         ; X  for _CDECL
      mov     edx,dword ptr[esp+8]         ; lX for _CDECL
      mov     ecx,dword ptr[esp+12]        ; Y  for _CDECL
      jmp     ?merge
_merge endp

public _lcmp                               ; _CDECL entry point

_lcmp proc
      mov     eax,dword ptr[esp+4]         ; pX
      mov     edx,dword ptr[esp+8]         ; pY
      jmp     ?lcmp
_lcmp endp

public ?merge

?merge PROC
      push    edi
      push    esi
      push    ebx                          ; for OPTLINK
      push    ebp
                                           ; set up limits
      sub     esp,8                        ; for 2 pointers
      mov     ebx,dword ptr[esp+40]        ; lY
      shl     ebx,2
      add     ebx,ecx
      mov     dword ptr[esp],ebx           ; end of sY

      mov     ebx,edx                      ; lX
      shl     ebx,2
      add     ebx,eax
      mov     dword ptr[esp+4],ebx         ; end of sX
      mov     edi,dword ptr[esp+44]        ; sZ

@LPHD:
      cmp     eax,dword ptr[esp+4]
      jae     @LAB5                        ; X exhausted
      cmp     ecx,dword ptr[esp]
      jae     @LAB4                        ; Y exhausted
      mov     ebx,dword ptr[eax]           ; pointer at x_i
      mov     ebp,dword ptr[ebx]           ; dereference, getting X.llcp
      mov     esi,dword ptr[ecx]           ; pointer at y_j
      cmp     ebp,dword ptr[esi]           ; deference, getting Y.llcp
      jg      @LAB1                        ; X.llcp > Y.llcp
      jl      @LAB2                        ; X.llcp < Y.llcp
      push    eax
      push    edx
      push    ecx                          ; save around function call
                                           ; eax, edx parameters to LCMP, OPTLINK linkage
      push    eax
      mov     eax, dword ptr [esp+68]      ; pointer to ncomp
      inc     dword ptr [eax]              ; update ncomp
      pop     eax
      mov     edx,ecx
      mov     ecx, dword ptr [esp+60]      ; lcol
      call    ?lcmp
      cmp     eax,0
      pop     ecx
      pop     edx
      pop     eax
      jg      @LAB2                        ; not in order, so output Y first
@LAB1:
      mov     ebx,dword ptr[eax]
      add     eax,4                        ; bump up pointer to sX
      jmp     @LAB3
@LAB2:
      mov     ebx,dword ptr[ecx]
      add     ecx,4                        ; bump up pointer to sY
@LAB3:
      mov     dword ptr[edi],ebx           ; set Z_k to A_i or B_j, as proper
      add     edi,4
      jmp     @LPHD                        ; loop
@LAB4:
      cmp     eax,dword ptr[esp+4]         ; X exhausted?
      je      @LAB5
      mov     ebx,dword ptr[eax]
      mov     dword ptr[edi],ebx           ; Z <- X
      add     eax,4
      add     edi,4
      jmp     @LAB4
@LAB5:
      add     esp,8
      pop     ebp
      pop     ebx
      pop     esi
      pop     edi
      ret
      ?merge  ENDP
;
      ; OPTLINK entry point
      ; i=lcmp(*sX,*sY,lcol). sX and sY are {int llcp; char *str;}
      ;
      ; compare strings starting at (sX.str)+sX.llcp, (sY.str)+sX.llcp
      ; stop at null or first non-matching character or at lcol.
      ; if first non-matching character of sX is greater, return 1, update sX.llcp
      ; Otherwise, return -1, update sY.llcp
      ; SHOULD NOT BE CALLED with sX.llcp different from sY.llcp
      ;

?lcmp PROC
      push    edi
      push    esi
      push    ebx                          ; for OPTLINK
      push    ebp                          ;
      mov     eax,dword ptr[eax]           ; dereference
      mov     edx,dword ptr[edx]
      mov     esi,dword ptr[eax]           ; sX.llcp
      mov     edi,dword ptr[edx]           ; sY.llcp
      mov     ebx,edx                      ; save offset of second par
      cmp     esi,edi
      jz      @BLB1
      xor     eax,eax
      jmp     @BLBL10
@BLB1:
      mov     ebp,dword  ptr[eax+4]        ; base of sX.str
      mov     edi,dword ptr [edx+4]        ; sY.str
      add     edi,esi                      ; first character of strY to compare
      lea     esi,dword ptr [ebp+esi]      ; sX.str+sX.llcp
      lea     ecx,dword ptr [ecx+ebp]      ; beyond lcol
      cld
@NEXT:
      cmp     esi,ecx
      jne     @NX                          ; at lcol+1 ?
      dec     esi
      jmp     @DONE
@NX:      
      mov     dl, byte ptr[esi]
      or      dl,dl
      jz      @DONE                        ; end of sx.str?
      cmpsb
      je      @NEXT
      dec     esi                          ; to compensate for autoincrement
@DONE:
      sub     esi,ebp                      ; new llcp of first string
      cmp     dl,[edi-1]
      jbe     @LESS                        ; first string is less
      mov     [eax],esi                    ; eax has not been changed
      mov     eax,1
      jmp     @BLBL10
@LESS:
      mov     [ebx],esi                    ; offset of second parameter
      mov     eax,-1
@BLBL10:
      pop     ebp
      pop     ebx
      pop     esi
      pop     edi
      ret
?lcmp endp
      END
