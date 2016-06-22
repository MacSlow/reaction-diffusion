global asmLaplaceSSE

section .data
	align 32
	matrixRow1: dd 0.05,   0.2,  0.05,  0.0
	matrixRow2: dd 0.2,   -1.0,  0.5,   0.0
	matrixRow3: dd 0.05,   0.2,  0.05,  0.0

; .05 * buf[offset-width-1]   .2 * buf[offset-width]   .05 * buf[offset-width+1]
; .2 * buf[offset-1]         -1. * buf[offset]         .2 * buf[offset+1]
; .05 * buf[offset+width-1]   .2 * buf[offset+width]   .05 * buf[offset+width+1]
;
;float laplaceA3 (BufferType& buffer, unsigned i, unsigned w) {
;    float sum = .0;
;    sum += get<0> (buffer[i]) * (-1);
;    sum += get<0> (buffer[i - 1]) * .2;
;    sum += get<0> (buffer[i + 1]) * .2;
;    sum += get<0> (buffer[i - w]) * .2;
;    sum += get<0> (buffer[i + w]) * .2;
;    sum += get<0> (buffer[i - 1 - w]) * .05;
;    sum += get<0> (buffer[i - 1 + w]) * .05;
;    sum += get<0> (buffer[i + 1 - w]) * .05;
;    sum += get<0> (buffer[i + 1 + w]) * .05;
;    return sum;
;}

section .text
; float asmLaplaceSSE (float* buffer,   // RDI
;                      unsigned offset, // offset into buffer in bytes in RSI
;                      unsigned width,  // width of row in bytes in RDX
;                      unsigned size)   // size of element in byte in RCX
; returns result in RAX/XMM0
asmLaplaceSSE:
	; prepare offsets
	mov r8, rdx
	imul r8, -1
	sub r8, 1
	mov r9, -1
	mov r10, rdx
	sub r10, 1

	imul r8, rcx
	imul r9, rcx
	imul r10, rcx

	add rdi, rsi

	; load 9 (12) floats from buffer
	movups xmm0, [rdi+r8]
	movups xmm1, [rdi+r9]
	movups xmm2, [rdi+r10]

	; load laplace-coefficients from memory
	movaps xmm3, [matrixRow1]
	movaps xmm4, [matrixRow2]
	movaps xmm5, [matrixRow3]

	; multiply 9 (12) floats
	mulps xmm0, xmm3
	mulps xmm1, xmm4
	mulps xmm2, xmm5

	; sum up all 9 (12) floats
	addps xmm0, xmm1
	addps xmm0, xmm2
	haddps xmm0, xmm0
	haddps xmm0, xmm0
	ret
