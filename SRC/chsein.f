      SUBROUTINE CHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, W, VL, LDVL, VR, LDVR, MM, M, WORK, RWORK, IFAILL, IFAILR, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EIGSRC, INITV, SIDE;
      int                INFO, LDH, LDVL, LDVR, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      int                IFAILL( * ), IFAILR( * );
      REAL               RWORK( * );
      COMPLEX            H( LDH, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      REAL               RZERO;
      const              RZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               BOTHV, FROMQR, LEFTV, NOINIT, RIGHTV;
      int                I, IINFO, K, KL, KLN, KR, KS, LDWORK;
      REAL               EPS3, HNORM, SMLNUM, ULP, UNFL;
      COMPLEX            CDUM, WK;
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      REAL               CLANHS, SLAMCH;
      // EXTERNAL LSAME, CLANHS, SLAMCH, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLAEIN, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) );
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters.

      BOTHV = LSAME( SIDE, 'B' );
      RIGHTV = LSAME( SIDE, 'R' ) || BOTHV;
      LEFTV = LSAME( SIDE, 'L' ) || BOTHV;

      FROMQR = LSAME( EIGSRC, 'Q' );

      NOINIT = LSAME( INITV, 'N' );

      // Set M to the number of columns required to store the selected
      // eigenvectors.

      M = 0;
      for (K = 1; K <= N; K++) { // 10
         IF( SELECT( K ) ) M = M + 1;
      } // 10

      INFO = 0;
      if ( !RIGHTV && !LEFTV ) {
         INFO = -1;
      } else if ( !FROMQR && !LSAME( EIGSRC, 'N' ) ) {
         INFO = -2;
      } else if ( !NOINIT && !LSAME( INITV, 'U' ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDH < MAX( 1, N ) ) {
         INFO = -7;
      } else if ( LDVL < 1 || ( LEFTV && LDVL < N ) ) {
         INFO = -10;
      } else if ( LDVR < 1 || ( RIGHTV && LDVR < N ) ) {
         INFO = -12;
      } else if ( MM < M ) {
         INFO = -13;
      }
      if ( INFO != 0 ) {
         xerbla('CHSEIN', -INFO );
         return;
      }

      // Quick return if possible.

      if (N == 0) RETURN;

      // Set machine-dependent constants.

      UNFL = SLAMCH( 'Safe minimum' );
      ULP = SLAMCH( 'Precision' );
      SMLNUM = UNFL*( N / ULP );

      LDWORK = N;

      KL = 1;
      KLN = 0;
      if ( FROMQR ) {
         KR = 0;
      } else {
         KR = N;
      }
      KS = 1;

      for (K = 1; K <= N; K++) { // 100
         if ( SELECT( K ) ) {

            // Compute eigenvector(s) corresponding to W(K).

            if ( FROMQR ) {

               // If affiliation of eigenvalues is known, check whether
               // the matrix splits.

               // Determine KL and KR such that 1 <= KL <= K <= KR <= N
               // and H(KL,KL-1) and H(KR+1,KR) are zero (or KL = 1 or
               // KR = N).

               // Then inverse iteration can be performed with the
               // submatrix H(KL:N,KL:N) for a left eigenvector, and with
               // the submatrix H(1:KR,1:KR) for a right eigenvector.

               DO 20 I = K, KL + 1, -1;
                  IF( H( I, I-1 ) == ZERO ) GO TO 30;
               } // 20
               } // 30
               KL = I;
               if ( K > KR ) {
                  for (I = K; I <= N - 1; I++) { // 40
                     IF( H( I+1, I ) == ZERO ) GO TO 50;
                  } // 40
                  } // 50
                  KR = I;
               }
            }

            if ( KL != KLN ) {
               KLN = KL;

               // Compute infinity-norm of submatrix H(KL:KR,KL:KR) if it
               // has not ben computed before.

               HNORM = CLANHS( 'I', KR-KL+1, H( KL, KL ), LDH, RWORK );
               if ( SISNAN( HNORM ) ) {
                  INFO = -6;
                  return;
               } else if ( (HNORM > RZERO) ) {
                  EPS3 = HNORM*ULP;
               } else {
                  EPS3 = SMLNUM;
               }
            }

            // Perturb eigenvalue if it is close to any previous
            // selected eigenvalues affiliated to the submatrix
            // H(KL:KR,KL:KR). Close roots are modified by EPS3.

            WK = W( K );
            } // 60
            DO 70 I = K - 1, KL, -1;
               if ( SELECT( I ) && CABS1( W( I )-WK ) < EPS3 ) {
                  WK = WK + EPS3;
                  GO TO 60;
               }
            } // 70
            W( K ) = WK;

            if ( LEFTV ) {

               // Compute left eigenvector.

               claein( false , NOINIT, N-KL+1, H( KL, KL ), LDH, WK, VL( KL, KS ), WORK, LDWORK, RWORK, EPS3, SMLNUM, IINFO );
               if ( IINFO > 0 ) {
                  INFO = INFO + 1;
                  IFAILL( KS ) = K;
               } else {
                  IFAILL( KS ) = 0;
               }
               for (I = 1; I <= KL - 1; I++) { // 80
                  VL( I, KS ) = ZERO;
               } // 80
            }
            if ( RIGHTV ) {

               // Compute right eigenvector.

               claein( true , NOINIT, KR, H, LDH, WK, VR( 1, KS ), WORK, LDWORK, RWORK, EPS3, SMLNUM, IINFO );
               if ( IINFO > 0 ) {
                  INFO = INFO + 1;
                  IFAILR( KS ) = K;
               } else {
                  IFAILR( KS ) = 0;
               }
               for (I = KR + 1; I <= N; I++) { // 90
                  VR( I, KS ) = ZERO;
               } // 90
            }
            KS = KS + 1;
         }
      } // 100

      return;

      // End of CHSEIN

      }
