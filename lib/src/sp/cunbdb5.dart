      void cunbdb5(final int M1, final int M2, final int N, final int X1, final int INCX1, final int X2, final int INCX2, final Matrix<double> Q1, final int LDQ1, final Matrix<double> Q2, final int LDQ2, final Array<double> WORK, final int LWORK, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2, N;
      Complex            Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*);
      // ..

      double               REALZERO;
      const              REALZERO = 0.0 ;
      Complex            ONE, ZERO;
      const              ONE = (1.0,0.0), ZERO = (0.0,0.0) ;
      int                CHILDINFO, I, J;
      double               EPS, NORM, SCL, SSQ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASSQ, CUNBDB6, CSCAL, XERBLA
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SCNRM2;
      // EXTERNAL SLAMCH, SCNRM2
      // ..
      // .. Intrinsic Function ..
      // INTRINSIC MAX

      // Test input arguments

      INFO = 0;
      if ( M1 < 0 ) {
         INFO = -1;
      } else if ( M2 < 0 ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( INCX1 < 1 ) {
         INFO = -5;
      } else if ( INCX2 < 1 ) {
         INFO = -7;
      } else if ( LDQ1 < max( 1, M1 ) ) {
         INFO = -9;
      } else if ( LDQ2 < max( 1, M2 ) ) {
         INFO = -11;
      } else if ( LWORK < N ) {
         INFO = -13;
      }

      if ( INFO != 0 ) {
         xerbla('CUNBDB5', -INFO );
         return;
      }

      EPS = SLAMCH( 'Precision' );

      // Project X onto the orthogonal complement of Q if X is nonzero

      SCL = REALZERO;
      SSQ = REALZERO;
      classq(M1, X1, INCX1, SCL, SSQ );
      classq(M2, X2, INCX2, SCL, SSQ );
      NORM = SCL * sqrt( SSQ );

      if ( NORM > N * EPS ) {
         // Scale vector to unit norm to avoid problems in the caller code.
         // Computing the reciprocal is undesirable but
         //  * xLASCL cannot be used because of the vector increments and
         //  * the round-off error has a negligible impact on
         //    orthogonalization.
         cscal(M1, ONE / NORM, X1, INCX1 );
         cscal(M2, ONE / NORM, X2, INCX2 );
         cunbdb6(M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK, CHILDINFO );

         // If the projection is nonzero, then return;

         if ( SCNRM2(M1,X1,INCX1) != REALZERO || SCNRM2(M2,X2,INCX2) != REALZERO ) {
            return;
         }
      }

      // Project each standard basis vector e_1,...,e_M1 in turn, stopping
      // when a nonzero projection is found

      for (I = 1; I <= M1; I++) {
         for (J = 1; J <= M1; J++) {
            X1[J] = ZERO;
         }
         X1[I] = ONE;
         for (J = 1; J <= M2; J++) {
            X2[J] = ZERO;
         }
         CALL CUNBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK, CHILDINFO )          IF( SCNRM2(M1,X1,INCX1) != REALZERO || SCNRM2(M2,X2,INCX2) != REALZERO ) THEN;
            return;
         }
      }

      // Project each standard basis vector e_(M1+1),...,e_(M1+M2) in turn,
      // stopping when a nonzero projection is found

      for (I = 1; I <= M2; I++) {
         for (J = 1; J <= M1; J++) {
            X1[J] = ZERO;
         }
         for (J = 1; J <= M2; J++) {
            X2[J] = ZERO;
         }
         X2[I] = ONE;
         CALL CUNBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, LDQ2, WORK, LWORK, CHILDINFO )          IF( SCNRM2(M1,X1,INCX1) != REALZERO || SCNRM2(M2,X2,INCX2) != REALZERO ) THEN;
            return;
         }
      }

      }
