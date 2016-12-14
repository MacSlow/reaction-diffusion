;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
;;
;; reaction-diffusion simulation in 2D
;;
;; Copyright 2016 Mirco Müller
;;
;; Author(s):
;;   Mirco "MacSlow" Müller <macslow@gmail.com>
;;
;; This program is free software: you can redistribute it and/or modify it
;; under the terms of the GNU General Public License version 3, as published
;; by the Free Software Foundation.
;;
;; This program is distributed in the hope that it will be useful, but
;; WITHOUT ANY WARRANTY; without even the implied warranties of
;; MERCHANTABILITY, SATISFACTORY QUALITY, or FITNESS FOR A PARTICULAR
;; PURPOSE.  See the GNU General Public License for more details.
;;
;; You should have received a copy of the GNU General Public License along
;; with this program.  If not, see <http:;;www.gnu.org/licenses/>.
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

global laplaceSSE
global calcA
global calcB

section .data
	align 16
	matrixRow1:	dd 0.05,   0.2,  0.05,  0.0
	matrixRow2:	dd 0.2,   -1.0,  0.2,   0.0
	matrixRow3:	dd 0.05,   0.2,  0.05,  0.0
	ones:		dd 1.0,    1.0,  1.0,   1.0

section .text
; float laplace (float* buffer,   // RDI
;                unsigned offset, // offset into buffer in bytes in RSI
;                unsigned width)  // width of row in bytes in RDX
; returns result in XMM0
laplaceSSE:
	; prepare offsets
	shl rsi, 2				; offset*4 because of sizeof(float) == 4
	add rdi, rsi			; rdi = buffer + offset
	sub rdi, 4				; rdi = buffer + offset - 4
	shl rdx, 2				; width*4 because of sizeof(float) == 4

	; load 9 (12) floats from buffer
	movups xmm1, [rdi]		; read four floats starting at buf+offset-1
	movups xmm2, [rdi+rdx]	; read four floats starting at buf+offset+width-1
	neg rdx
	movups xmm0, [rdi+rdx]	; read four floats starting at buf+offset-width-1

	; load laplace-coefficients from memory
	movaps xmm3, [rel matrixRow1]
	movaps xmm4, [rel matrixRow2]
	movaps xmm5, [rel matrixRow3]

	; multiply 9 (12) floats
	mulps xmm0, xmm3
	mulps xmm1, xmm4
	mulps xmm2, xmm5

	; sum up all 9 (12) floats
	addps xmm0, xmm1
	addps xmm0, xmm2
	haddps xmm0, xmm0		; horizontal add is cool, but costly
	haddps xmm0, xmm0		; horizontal add is cool, but costly

	ret

; FIXME: still unfinished
calcA:
;void calcA (unsigned index,  // RDI
;			 unsigned width,  // RSI
;			 float* aScratch, // RDX  ... xmm7
;			 float* aBuffer,  // RCX  ... xmm2
;			 float* bBuffer,  // R8   ... xmm4
;			 float* feed,     // R9   ... xmm5
;			 float dA,        //          xmm0
;			 float dt)        //          xmm1
; returns result into aScratch
; aScratch = a    + (dA   * laplaceA - a*b*b + feed  * (1 - a)) * dt
; xmm7     = xmm2 + (xmm0 * xmm3     - xmm4  + xmm5  *  xmm6)   * xmm1
; fma1:              xmm0 * xmm3     - xmm4
; fma2:              fma1                    + xmm5  *  xmm6
; fma3:      xmm2 +  fma2                                       * xmm1
;
; laplace:
;    rdi - buffer
;    rsi - offset/index
;    rdx - width
;
; local to laplace:
;    r8, r9, xmm0, xmm1, xmm2, xmm3, xmm4, xmm5
;
	push rbp
	mov rbp, rsp

	shl rdi, 2
	shl rsi, 2
	add rdi, rsi
	movups xmm7, [rdx]		; aScratch
	movups xmm2, [rcx]		; a
	movups xmm4, [r8]		; b
	movups xmm5, [r9]		; feed

	movaps xmm6, [rel ones]		; 1 - a
	subps xmm6, xmm2

	mulps xmm4, xmm4		; a * b * b
	mulps xmm4, xmm2

	mov rsp, rbp
	pop rbp
	ret

; FIXME: still unfinished
calcB:
;void calcB (unsigned index,  // RDI
;			 unsigned width,  // RSI
;			 float* bScratch, // RDX         ... xmm7
;			 float* aBuffer,  // RCX         ... xmm4
;			 float* bBuffer,  // R8          ... xmm2
;			 float* feed,     // R9          ... xmm5
;		 	 float* kill,     // stack argl! ... xmm6
;			 float dB,        //                 xmm0
;			 float dt)        //                 xmm1
; returns result into bScratch
; bScratch = b    + (dB   * laplaceB + a*b*b - (feed + kill) *    b) * dt
; xmm7     = xmm2 + (xmm0 * xmm3     + xmm4  - (xmm5 + xmm6) * xmm2) * xmm1
;                                      xmm9
; xmm0:              xmm0 * xmm3     + xmm9
; ....:                                         xmm8
; xmm0:              xmm0                    -  xmm8         * xmm2
; fma3:      xmm2 +  xmm0                                            * xmm1

	push rbp
	mov rbp, rsp

	mov rax, rdi			; save local arguments
	mov r10, rsi
	mov r11, rdx
	mov rdx, rsi			; prepare arguments for laplace ()
	mov rsi, rdi
	mov rdi, r8
	call laplaceSSE
	movaps xmm3, xmm0 		; returns result in xmm3
	mov rdi, rax			; restore local arguments
	mov rsi, r10
	mov rdx, r11

	mov rax, [rbp+16]		; grab pointer to kill from stack

	movaps xmm4, [rcx]		; a
	movaps xmm2, [r8]		; b
	movaps xmm5, [r9]		; feed
	movaps xmm6, [rax]		; kill

	movaps xmm9, xmm2		; a * b * b
	mulps xmm9, xmm9
	mulps xmm9, xmm4

	movaps xmm8, xmm5		; feed + kill
	addps xmm8, xmm6

	mulps xmm0, xmm3
	addps xmm0, xmm9
	mulps xmm8, xmm2
	subps xmm0, xmm8
	mulps xmm0, xmm1
	addps xmm0, xmm2

	movups [rdx], xmm0

	mov rsp, rbp
	pop rbp
	ret
