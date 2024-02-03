      SUBROUTINE SHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, WR, WI, VL, LDVL, VR, LDVR, MM, M, WORK, IFAILL, IFAILR, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EIGSRC, INITV, SIDE;
      int                INFO, LDH, LDVL, LDVR, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      int                IFAILL( * ), IFAILR( * );
      REAL               H( LDH, * ), VL( LDVL, * ), VR( LDVR, * ), WI( * ), WORK( * ), WR( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               BOTHV, FROMQR, LEFTV, NOINIT, PAIR, RIGHTV;
      int                I, IINFO, K, KL, KLN, KR, KSI, KSR, LDWORK;
      REAL               BIGNUM, EPS3, HNORM, SMLNUM, ULP, UNFL, WKI, WKR
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      REAL               SLAMCH, SLANHS
      // EXTERNAL LSAME, SLAMCH, SLANHS, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAEIN, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters.

      BOTHV = LSAME( SIDE, 'B' )
      RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
      LEFTV = LSAME( SIDE, 'L' ) .OR. BOTHV

      FROMQR = LSAME( EIGSRC, 'Q' )

      NOINIT = LSAME( INITV, 'N' )

      // Set M to the number of columns required to store the selected
      // eigenvectors, and standardize the array SELECT.

      M = 0
      PAIR = false;
      for (K = 1; K <= N; K++) { // 10
         if ( PAIR ) {
            PAIR = false;
            SELECT( K ) = false;
         } else {
            if ( WI( K ).EQ.ZERO ) {
               IF( SELECT( K ) ) M = M + 1
            } else {
               PAIR = true;
               if ( SELECT( K ) .OR. SELECT( K+1 ) ) {
                  SELECT( K ) = true;
                  M = M + 2
               }
            }
         }
      } // 10

      INFO = 0
      if ( .NOT.RIGHTV .AND. .NOT.LEFTV ) {
         INFO = -1
      } else if ( .NOT.FROMQR .AND. .NOT.LSAME( EIGSRC, 'N' ) ) {
         INFO = -2
      } else if ( .NOT.NOINIT .AND. .NOT.LSAME( INITV, 'U' ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -5
      } else if ( LDH.LT.MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDVL.LT.1 .OR. ( LEFTV .AND. LDVL.LT.N ) ) {
         INFO = -11
      } else if ( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) {
         INFO = -13
      } else if ( MM.LT.M ) {
         INFO = -14
      }
      if ( INFO.NE.0 ) {
         xerbla('SHSEIN', -INFO );
         RETURN
      }

      // Quick return if possible.

      if (N.EQ.0) RETURN;

      // Set machine-dependent constants.

      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )
      BIGNUM = ( ONE-ULP ) / SMLNUM

      LDWORK = N + 1

      KL = 1
      KLN = 0
      if ( FROMQR ) {
         KR = 0
      } else {
         KR = N
      }
      KSR = 1

      for (K = 1; K <= N; K++) { // 120
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

               DO 20 I = K, KL + 1, -1
                  IF( H( I, I-1 ).EQ.ZERO ) GO TO 30
               } // 20
               } // 30
               KL = I
               if ( K.GT.KR ) {
                  for (I = K; I <= N - 1; I++) { // 40
                     IF( H( I+1, I ).EQ.ZERO ) GO TO 50
                  } // 40
                  } // 50
                  KR = I
               }
            }

            if ( KL.NE.KLN ) {
               KLN = KL

               // Compute infinity-norm of submatrix H(KL:KR,KL:KR) if it
               // has not ben computed before.

               HNORM = SLANHS( 'I', KR-KL+1, H( KL, KL ), LDH, WORK )
               if ( SISNAN( HNORM ) ) {
                  INFO = -6
                  RETURN
               } else if ( HNORM.GT.ZERO ) {
                  EPS3 = HNORM*ULP
               } else {
                  EPS3 = SMLNUM
               }
            }

            // Perturb eigenvalue if it is close to any previous
            // selected eigenvalues affiliated to the submatrix
            // H(KL:KR,KL:KR). Close roots are modified by EPS3.

            WKR = WR( K )
            WKI = WI( K )
            } // 60
            DO 70 I = K - 1, KL, -1
               if ( SELECT( I ) .AND. ABS( WR( I )-WKR )+ ABS( WI( I )-WKI ).LT.EPS3 ) {
                  WKR = WKR + EPS3
                  GO TO 60
               }
            } // 70
            WR( K ) = WKR

            PAIR = WKI.NE.ZERO
            if ( PAIR ) {
               KSI = KSR + 1
            } else {
               KSI = KSR
            }
            if ( LEFTV ) {

               // Compute left eigenvector.

               slaein( false , NOINIT, N-KL+1, H( KL, KL ), LDH, WKR, WKI, VL( KL, KSR ), VL( KL, KSI ), WORK, LDWORK, WORK( N*N+N+1 ), EPS3, SMLNUM, BIGNUM, IINFO );
               if ( IINFO.GT.0 ) {
                  if ( PAIR ) {
                     INFO = INFO + 2
                  } else {
                     INFO = INFO + 1
                  }
                  IFAILL( KSR ) = K
                  IFAILL( KSI ) = K
               } else {
                  IFAILL( KSR ) = 0
                  IFAILL( KSI ) = 0
               }
               for (I = 1; I <= KL - 1; I++) { // 80
                  VL( I, KSR ) = ZERO
               } // 80
               if ( PAIR ) {
                  for (I = 1; I <= KL - 1; I++) { // 90
                     VL( I, KSI ) = ZERO
                  } // 90
               }
            }
            if ( RIGHTV ) {

               // Compute right eigenvector.

               slaein( true , NOINIT, KR, H, LDH, WKR, WKI, VR( 1, KSR ), VR( 1, KSI ), WORK, LDWORK, WORK( N*N+N+1 ), EPS3, SMLNUM, BIGNUM, IINFO );
               if ( IINFO.GT.0 ) {
                  if ( PAIR ) {
                     INFO = INFO + 2
                  } else {
                     INFO = INFO + 1
                  }
                  IFAILR( KSR ) = K
                  IFAILR( KSI ) = K
               } else {
                  IFAILR( KSR ) = 0
                  IFAILR( KSI ) = 0
               }
               for (I = KR + 1; I <= N; I++) { // 100
                  VR( I, KSR ) = ZERO
               } // 100
               if ( PAIR ) {
                  for (I = KR + 1; I <= N; I++) { // 110
                     VR( I, KSI ) = ZERO
                  } // 110
               }
            }

            if ( PAIR ) {
               KSR = KSR + 2
            } else {
               KSR = KSR + 1
            }
         }
      } // 120

      RETURN

      // End of SHSEIN

      }
