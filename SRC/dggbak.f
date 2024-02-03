      SUBROUTINE DGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V, LDV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOB, SIDE;
      int                IHI, ILO, INFO, LDV, M, N;
      // ..
      // .. Array Arguments ..
      double             LSCALE( * ), RSCALE( * ), V( LDV, * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LEFTV, RIGHTV;
      int                I, K;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, INT
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      RIGHTV = LSAME( SIDE, 'R' )
      LEFTV = LSAME( SIDE, 'L' )

      INFO = 0
      if ( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) {
         INFO = -1
      } else if ( .NOT.RIGHTV .AND. .NOT.LEFTV ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( ILO.LT.1 ) {
         INFO = -4
      } else if ( N.EQ.0 .AND. IHI.EQ.0 .AND. ILO.NE.1 ) {
         INFO = -4
      } else if ( N.GT.0 .AND. ( IHI.LT.ILO .OR. IHI.GT.MAX( 1, N ) ) ) {
         INFO = -5
      } else if ( N.EQ.0 .AND. ILO.EQ.1 .AND. IHI.NE.0 ) {
         INFO = -5
      } else if ( M.LT.0 ) {
         INFO = -8
      } else if ( LDV.LT.MAX( 1, N ) ) {
         INFO = -10
      }
      if ( INFO.NE.0 ) {
         xerbla('DGGBAK', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN       IF( M.EQ.0 ) RETURN       IF( LSAME( JOB, 'N' ) ) RETURN

      IF( ILO.EQ.IHI ) GO TO 30

      // Backward balance

      if ( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) {

         // Backward transformation on right eigenvectors

         if ( RIGHTV ) {
            DO 10 I = ILO, IHI
               dscal(M, RSCALE( I ), V( I, 1 ), LDV );
   10       CONTINUE
         }

         // Backward transformation on left eigenvectors

         if ( LEFTV ) {
            DO 20 I = ILO, IHI
               dscal(M, LSCALE( I ), V( I, 1 ), LDV );
   20       CONTINUE
         }
      }

      // Backward permutation

   30 CONTINUE
      if ( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) {

         // Backward permutation on right eigenvectors

         if ( RIGHTV ) {
            IF( ILO.EQ.1 ) GO TO 50

            DO 40 I = ILO - 1, 1, -1
               K = INT(RSCALE( I ))
               IF( K.EQ.I ) GO TO 40
               dswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
   40       CONTINUE

   50       CONTINUE
            IF( IHI.EQ.N ) GO TO 70
            DO 60 I = IHI + 1, N
               K = INT(RSCALE( I ))
               IF( K.EQ.I ) GO TO 60
               dswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
   60       CONTINUE
         }

         // Backward permutation on left eigenvectors

   70    CONTINUE
         if ( LEFTV ) {
            IF( ILO.EQ.1 ) GO TO 90
            DO 80 I = ILO - 1, 1, -1
               K = INT(LSCALE( I ))
               IF( K.EQ.I ) GO TO 80
               dswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
   80       CONTINUE

   90       CONTINUE
            IF( IHI.EQ.N ) GO TO 110
            DO 100 I = IHI + 1, N
               K = INT(LSCALE( I ))
               IF( K.EQ.I ) GO TO 100
               dswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
  100       CONTINUE
         }
      }

  110 CONTINUE

      RETURN

      // End of DGGBAK

      }
