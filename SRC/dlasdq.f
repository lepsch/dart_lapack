      SUBROUTINE DLASDQ( UPLO, SQRE, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, WORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU, SQRE;
      // ..
      // .. Array Arguments ..
      double             C( LDC, * ), D( * ), E( * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               ROTATE;
      int                I, ISUB, IUPLO, J, NP1, SQRE1;
      double             CS, R, SMIN, SN;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DBDSQR, DLARTG, DLASR, DSWAP, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      IUPLO = 0
      IF( LSAME( UPLO, 'U' ) ) IUPLO = 1       IF( LSAME( UPLO, 'L' ) ) IUPLO = 2
      if ( IUPLO.EQ.0 ) {
         INFO = -1
      } else if ( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( NCVT.LT.0 ) {
         INFO = -4
      } else if ( NRU.LT.0 ) {
         INFO = -5
      } else if ( NCC.LT.0 ) {
         INFO = -6
      } else if ( ( NCVT.EQ.0 .AND. LDVT.LT.1 ) .OR. ( NCVT.GT.0 .AND. LDVT.LT.MAX( 1, N ) ) ) {
         INFO = -10
      } else if ( LDU.LT.MAX( 1, NRU ) ) {
         INFO = -12
      } else if ( ( NCC.EQ.0 .AND. LDC.LT.1 ) .OR. ( NCC.GT.0 .AND. LDC.LT.MAX( 1, N ) ) ) {
         INFO = -14
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DLASDQ', -INFO )
         RETURN
      }
      IF( N.EQ.0 ) RETURN

      // ROTATE is true if any singular vectors desired, false otherwise

      ROTATE = ( NCVT.GT.0 ) .OR. ( NRU.GT.0 ) .OR. ( NCC.GT.0 )
      NP1 = N + 1
      SQRE1 = SQRE

      // If matrix non-square upper bidiagonal, rotate to be lower
      // bidiagonal.  The rotations are on the right.

      if ( ( IUPLO.EQ.1 ) .AND. ( SQRE1.EQ.1 ) ) {
         DO 10 I = 1, N - 1
            CALL DLARTG( D( I ), E( I ), CS, SN, R )
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            if ( ROTATE ) {
               WORK( I ) = CS
               WORK( N+I ) = SN
            }
   10    CONTINUE
         CALL DLARTG( D( N ), E( N ), CS, SN, R )
         D( N ) = R
         E( N ) = ZERO
         if ( ROTATE ) {
            WORK( N ) = CS
            WORK( N+N ) = SN
         }
         IUPLO = 2
         SQRE1 = 0

         // Update singular vectors if desired.

         IF( NCVT.GT.0 ) CALL DLASR( 'L', 'V', 'F', NP1, NCVT, WORK( 1 ), WORK( NP1 ), VT, LDVT )
      }

      // If matrix lower bidiagonal, rotate to be upper bidiagonal
      // by applying Givens rotations on the left.

      if ( IUPLO.EQ.2 ) {
         DO 20 I = 1, N - 1
            CALL DLARTG( D( I ), E( I ), CS, SN, R )
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            if ( ROTATE ) {
               WORK( I ) = CS
               WORK( N+I ) = SN
            }
   20    CONTINUE

         // If matrix (N+1)-by-N lower bidiagonal, one additional
         // rotation is needed.

         if ( SQRE1.EQ.1 ) {
            CALL DLARTG( D( N ), E( N ), CS, SN, R )
            D( N ) = R
            if ( ROTATE ) {
               WORK( N ) = CS
               WORK( N+N ) = SN
            }
         }

         // Update singular vectors if desired.

         if ( NRU.GT.0 ) {
            if ( SQRE1.EQ.0 ) {
               CALL DLASR( 'R', 'V', 'F', NRU, N, WORK( 1 ), WORK( NP1 ), U, LDU )
            } else {
               CALL DLASR( 'R', 'V', 'F', NRU, NP1, WORK( 1 ), WORK( NP1 ), U, LDU )
            }
         }
         if ( NCC.GT.0 ) {
            if ( SQRE1.EQ.0 ) {
               CALL DLASR( 'L', 'V', 'F', N, NCC, WORK( 1 ), WORK( NP1 ), C, LDC )
            } else {
               CALL DLASR( 'L', 'V', 'F', NP1, NCC, WORK( 1 ), WORK( NP1 ), C, LDC )
            }
         }
      }

      // Call DBDSQR to compute the SVD of the reduced real
      // N-by-N upper bidiagonal matrix.

      CALL DBDSQR( 'U', N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, WORK, INFO )

      // Sort the singular values into ascending order (insertion sort on
      // singular values, but only one transposition per singular vector)

      DO 40 I = 1, N

         // Scan for smallest D(I).

         ISUB = I
         SMIN = D( I )
         DO 30 J = I + 1, N
            if ( D( J ).LT.SMIN ) {
               ISUB = J
               SMIN = D( J )
            }
   30    CONTINUE
         if ( ISUB.NE.I ) {

            // Swap singular values and vectors.

            D( ISUB ) = D( I )
            D( I ) = SMIN
            IF( NCVT.GT.0 ) CALL DSWAP( NCVT, VT( ISUB, 1 ), LDVT, VT( I, 1 ), LDVT )             IF( NRU.GT.0 ) CALL DSWAP( NRU, U( 1, ISUB ), 1, U( 1, I ), 1 )             IF( NCC.GT.0 ) CALL DSWAP( NCC, C( ISUB, 1 ), LDC, C( I, 1 ), LDC )
         }
   40 CONTINUE

      RETURN

      // End of DLASDQ

      }
