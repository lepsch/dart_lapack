      SUBROUTINE SGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V, LDV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOB, SIDE;
      int                IHI, ILO, INFO, LDV, M, N;
      // ..
      // .. Array Arguments ..
      REAL               LSCALE( * ), RSCALE( * ), V( LDV, * )
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
      // EXTERNAL SSCAL, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
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
         xerbla('SGGBAK', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N.EQ.0) RETURN       IF( M.EQ.0 ) RETURN       IF( LSAME( JOB, 'N' ) ) RETURN;

      if (ILO.EQ.IHI) GO TO 30;

      // Backward balance

      if ( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) {

         // Backward transformation on right eigenvectors

         if ( RIGHTV ) {
            for (I = ILO; I <= IHI; I++) { // 10
               sscal(M, RSCALE( I ), V( I, 1 ), LDV );
            } // 10
         }

         // Backward transformation on left eigenvectors

         if ( LEFTV ) {
            for (I = ILO; I <= IHI; I++) { // 20
               sscal(M, LSCALE( I ), V( I, 1 ), LDV );
            } // 20
         }
      }

      // Backward permutation

      } // 30
      if ( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) {

         // Backward permutation on right eigenvectors

         if ( RIGHTV ) {
            if (ILO.EQ.1) GO TO 50;

            DO 40 I = ILO - 1, 1, -1
               K = INT( RSCALE( I ) )
               if (K.EQ.I) GO TO 40;
               sswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
            } // 40

            } // 50
            if (IHI.EQ.N) GO TO 70;
            for (I = IHI + 1; I <= N; I++) { // 60
               K = INT( RSCALE( I ) )
               if (K.EQ.I) GO TO 60;
               sswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
            } // 60
         }

         // Backward permutation on left eigenvectors

         } // 70
         if ( LEFTV ) {
            if (ILO.EQ.1) GO TO 90;
            DO 80 I = ILO - 1, 1, -1
               K = INT( LSCALE( I ) )
               if (K.EQ.I) GO TO 80;
               sswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
            } // 80

            } // 90
            if (IHI.EQ.N) GO TO 110;
            for (I = IHI + 1; I <= N; I++) { // 100
               K = INT( LSCALE( I ) )
               if (K.EQ.I) GO TO 100;
               sswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
            } // 100
         }
      }

      } // 110

      RETURN

      // End of SGGBAK

      }
