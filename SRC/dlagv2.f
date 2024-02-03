      SUBROUTINE DLAGV2( A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, CSL, SNL, CSR, SNR )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LDB
      DOUBLE PRECISION   CSL, CSR, SNL, SNR
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( 2 ), ALPHAR( 2 ), B( LDB, * ), BETA( 2 )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   ANORM, ASCALE, BNORM, BSCALE, H1, H2, H3, QQ, R, RR, SAFMIN, SCALE1, SCALE2, T, ULP, WI, WR1, WR2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAG2, DLARTG, DLASV2, DROT
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           DLAMCH, DLAPY2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
      SAFMIN = DLAMCH( 'S' )
      ULP = DLAMCH( 'P' )
*
*     Scale A
*
      ANORM = MAX( ABS( A( 1, 1 ) )+ABS( A( 2, 1 ) ), ABS( A( 1, 2 ) )+ABS( A( 2, 2 ) ), SAFMIN )
      ASCALE = ONE / ANORM
      A( 1, 1 ) = ASCALE*A( 1, 1 )
      A( 1, 2 ) = ASCALE*A( 1, 2 )
      A( 2, 1 ) = ASCALE*A( 2, 1 )
      A( 2, 2 ) = ASCALE*A( 2, 2 )
*
*     Scale B
*
      BNORM = MAX( ABS( B( 1, 1 ) ), ABS( B( 1, 2 ) )+ABS( B( 2, 2 ) ), SAFMIN )
      BSCALE = ONE / BNORM
      B( 1, 1 ) = BSCALE*B( 1, 1 )
      B( 1, 2 ) = BSCALE*B( 1, 2 )
      B( 2, 2 ) = BSCALE*B( 2, 2 )
*
*     Check if A can be deflated
*
      IF( ABS( A( 2, 1 ) ).LE.ULP ) THEN
         CSL = ONE
         SNL = ZERO
         CSR = ONE
         SNR = ZERO
         A( 2, 1 ) = ZERO
         B( 2, 1 ) = ZERO
         WI = ZERO
*
*     Check if B is singular
*
      ELSE IF( ABS( B( 1, 1 ) ).LE.ULP ) THEN
         CALL DLARTG( A( 1, 1 ), A( 2, 1 ), CSL, SNL, R )
         CSR = ONE
         SNR = ZERO
         CALL DROT( 2, A( 1, 1 ), LDA, A( 2, 1 ), LDA, CSL, SNL )
         CALL DROT( 2, B( 1, 1 ), LDB, B( 2, 1 ), LDB, CSL, SNL )
         A( 2, 1 ) = ZERO
         B( 1, 1 ) = ZERO
         B( 2, 1 ) = ZERO
         WI = ZERO
*
      ELSE IF( ABS( B( 2, 2 ) ).LE.ULP ) THEN
         CALL DLARTG( A( 2, 2 ), A( 2, 1 ), CSR, SNR, T )
         SNR = -SNR
         CALL DROT( 2, A( 1, 1 ), 1, A( 1, 2 ), 1, CSR, SNR )
         CALL DROT( 2, B( 1, 1 ), 1, B( 1, 2 ), 1, CSR, SNR )
         CSL = ONE
         SNL = ZERO
         A( 2, 1 ) = ZERO
         B( 2, 1 ) = ZERO
         B( 2, 2 ) = ZERO
         WI = ZERO
*
      ELSE
*
*        B is nonsingular, first compute the eigenvalues of (A,B)
*
         CALL DLAG2( A, LDA, B, LDB, SAFMIN, SCALE1, SCALE2, WR1, WR2, WI )
*
         IF( WI.EQ.ZERO ) THEN
*
*           two real eigenvalues, compute s*A-w*B
*
            H1 = SCALE1*A( 1, 1 ) - WR1*B( 1, 1 )
            H2 = SCALE1*A( 1, 2 ) - WR1*B( 1, 2 )
            H3 = SCALE1*A( 2, 2 ) - WR1*B( 2, 2 )
*
            RR = DLAPY2( H1, H2 )
            QQ = DLAPY2( SCALE1*A( 2, 1 ), H3 )
*
            IF( RR.GT.QQ ) THEN
*
*              find right rotation matrix to zero 1,1 element of
*              (sA - wB)
*
               CALL DLARTG( H2, H1, CSR, SNR, T )
*
            ELSE
*
*              find right rotation matrix to zero 2,1 element of
*              (sA - wB)
*
               CALL DLARTG( H3, SCALE1*A( 2, 1 ), CSR, SNR, T )
*
            END IF
*
            SNR = -SNR
            CALL DROT( 2, A( 1, 1 ), 1, A( 1, 2 ), 1, CSR, SNR )
            CALL DROT( 2, B( 1, 1 ), 1, B( 1, 2 ), 1, CSR, SNR )
*
*           compute inf norms of A and B
*
            H1 = MAX( ABS( A( 1, 1 ) )+ABS( A( 1, 2 ) ), ABS( A( 2, 1 ) )+ABS( A( 2, 2 ) ) )             H2 = MAX( ABS( B( 1, 1 ) )+ABS( B( 1, 2 ) ), ABS( B( 2, 1 ) )+ABS( B( 2, 2 ) ) )
*
            IF( ( SCALE1*H1 ).GE.ABS( WR1 )*H2 ) THEN
*
*              find left rotation matrix Q to zero out B(2,1)
*
               CALL DLARTG( B( 1, 1 ), B( 2, 1 ), CSL, SNL, R )
*
            ELSE
*
*              find left rotation matrix Q to zero out A(2,1)
*
               CALL DLARTG( A( 1, 1 ), A( 2, 1 ), CSL, SNL, R )
*
            END IF
*
            CALL DROT( 2, A( 1, 1 ), LDA, A( 2, 1 ), LDA, CSL, SNL )
            CALL DROT( 2, B( 1, 1 ), LDB, B( 2, 1 ), LDB, CSL, SNL )
*
            A( 2, 1 ) = ZERO
            B( 2, 1 ) = ZERO
*
         ELSE
*
*           a pair of complex conjugate eigenvalues
*           first compute the SVD of the matrix B
*
            CALL DLASV2( B( 1, 1 ), B( 1, 2 ), B( 2, 2 ), R, T, SNR, CSR, SNL, CSL )
*
*           Form (A,B) := Q(A,B)Z**T where Q is left rotation matrix and
*           Z is right rotation matrix computed from DLASV2
*
            CALL DROT( 2, A( 1, 1 ), LDA, A( 2, 1 ), LDA, CSL, SNL )
            CALL DROT( 2, B( 1, 1 ), LDB, B( 2, 1 ), LDB, CSL, SNL )
            CALL DROT( 2, A( 1, 1 ), 1, A( 1, 2 ), 1, CSR, SNR )
            CALL DROT( 2, B( 1, 1 ), 1, B( 1, 2 ), 1, CSR, SNR )
*
            B( 2, 1 ) = ZERO
            B( 1, 2 ) = ZERO
*
         END IF
*
      END IF
*
*     Unscaling
*
      A( 1, 1 ) = ANORM*A( 1, 1 )
      A( 2, 1 ) = ANORM*A( 2, 1 )
      A( 1, 2 ) = ANORM*A( 1, 2 )
      A( 2, 2 ) = ANORM*A( 2, 2 )
      B( 1, 1 ) = BNORM*B( 1, 1 )
      B( 2, 1 ) = BNORM*B( 2, 1 )
      B( 1, 2 ) = BNORM*B( 1, 2 )
      B( 2, 2 ) = BNORM*B( 2, 2 )
*
      IF( WI.EQ.ZERO ) THEN
         ALPHAR( 1 ) = A( 1, 1 )
         ALPHAR( 2 ) = A( 2, 2 )
         ALPHAI( 1 ) = ZERO
         ALPHAI( 2 ) = ZERO
         BETA( 1 ) = B( 1, 1 )
         BETA( 2 ) = B( 2, 2 )
      ELSE
         ALPHAR( 1 ) = ANORM*WR1 / SCALE1 / BNORM
         ALPHAI( 1 ) = ANORM*WI / SCALE1 / BNORM
         ALPHAR( 2 ) = ALPHAR( 1 )
         ALPHAI( 2 ) = -ALPHAI( 1 )
         BETA( 1 ) = ONE
         BETA( 2 ) = ONE
      END IF
*
      RETURN
*
*     End of DLAGV2
*
      END
