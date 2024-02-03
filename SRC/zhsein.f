      SUBROUTINE ZHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, W, VL, LDVL, VR, LDVR, MM, M, WORK, RWORK, IFAILL, IFAILR, INFO )

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
      double             RWORK( * );
      COMPLEX*16         H( LDH, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO
      const              ZERO = ( 0.0D+0, 0.0D+0 ) ;
      double             RZERO;
      const              RZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               BOTHV, FROMQR, LEFTV, NOINIT, RIGHTV;
      int                I, IINFO, K, KL, KLN, KR, KS, LDWORK;
      double             EPS3, HNORM, SMLNUM, ULP, UNFL;
      COMPLEX*16         CDUM, WK
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      double             DLAMCH, ZLANHS;
      // EXTERNAL LSAME, DLAMCH, ZLANHS, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLAEIN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters.

      BOTHV = LSAME( SIDE, 'B' )
      RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
      LEFTV = LSAME( SIDE, 'L' ) .OR. BOTHV

      FROMQR = LSAME( EIGSRC, 'Q' )

      NOINIT = LSAME( INITV, 'N' )

      // Set M to the number of columns required to store the selected
      // eigenvectors.

      M = 0
      for (K = 1; K <= N; K++) { // 10
         IF( SELECT( K ) ) M = M + 1
   10 CONTINUE

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
         INFO = -10
      } else if ( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) {
         INFO = -12
      } else if ( MM.LT.M ) {
         INFO = -13
      }
      if ( INFO.NE.0 ) {
         xerbla('ZHSEIN', -INFO );
         RETURN
      }

      // Quick return if possible.

      IF( N.EQ.0 ) RETURN

      // Set machine-dependent constants.

      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )

      LDWORK = N

      KL = 1
      KLN = 0
      if ( FROMQR ) {
         KR = 0
      } else {
         KR = N
      }
      KS = 1

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

               DO 20 I = K, KL + 1, -1
                  IF( H( I, I-1 ).EQ.ZERO ) GO TO 30
   20          CONTINUE
   30          CONTINUE
               KL = I
               if ( K.GT.KR ) {
                  DO 40 I = K, N - 1
                     IF( H( I+1, I ).EQ.ZERO ) GO TO 50
   40             CONTINUE
   50             CONTINUE
                  KR = I
               }
            }

            if ( KL.NE.KLN ) {
               KLN = KL

               // Compute infinity-norm of submatrix H(KL:KR,KL:KR) if it
               // has not ben computed before.

               HNORM = ZLANHS( 'I', KR-KL+1, H( KL, KL ), LDH, RWORK )
               if ( DISNAN( HNORM ) ) {
                  INFO = -6
                  RETURN
               } else if ( HNORM.GT.RZERO ) {
                  EPS3 = HNORM*ULP
               } else {
                  EPS3 = SMLNUM
               }
            }

            // Perturb eigenvalue if it is close to any previous
            // selected eigenvalues affiliated to the submatrix
            // H(KL:KR,KL:KR). Close roots are modified by EPS3.

            WK = W( K )
   60       CONTINUE
            DO 70 I = K - 1, KL, -1
               if ( SELECT( I ) .AND. CABS1( W( I )-WK ).LT.EPS3 ) {
                  WK = WK + EPS3
                  GO TO 60
               }
   70       CONTINUE
            W( K ) = WK

            if ( LEFTV ) {

               // Compute left eigenvector.

               zlaein(.FALSE., NOINIT, N-KL+1, H( KL, KL ), LDH, WK, VL( KL, KS ), WORK, LDWORK, RWORK, EPS3, SMLNUM, IINFO );
               if ( IINFO.GT.0 ) {
                  INFO = INFO + 1
                  IFAILL( KS ) = K
               } else {
                  IFAILL( KS ) = 0
               }
               DO 80 I = 1, KL - 1
                  VL( I, KS ) = ZERO
   80          CONTINUE
            }
            if ( RIGHTV ) {

               // Compute right eigenvector.

               zlaein(.TRUE., NOINIT, KR, H, LDH, WK, VR( 1, KS ), WORK, LDWORK, RWORK, EPS3, SMLNUM, IINFO );
               if ( IINFO.GT.0 ) {
                  INFO = INFO + 1
                  IFAILR( KS ) = K
               } else {
                  IFAILR( KS ) = 0
               }
               DO 90 I = KR + 1, N
                  VR( I, KS ) = ZERO
   90          CONTINUE
            }
            KS = KS + 1
         }
  100 CONTINUE

      RETURN

      // End of ZHSEIN

      }
