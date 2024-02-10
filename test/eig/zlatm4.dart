      void zlatm4(ITYPE, N, NZ1, NZ2, RSIGN, AMAGN, RCOND, TRIANG, IDIST, ISEED, A, LDA ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               RSIGN;
      int                IDIST, ITYPE, LDA, N, NZ1, NZ2;
      double             AMAGN, RCOND, TRIANG;
      int                ISEED( 4 );
      Complex         A( LDA, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                I, ISDB, ISDE, JC, JD, JR, K, KBEG, KEND, KLEN;
      double             ALPHA;
      Complex         CTEMP;
      // ..
      // .. External Functions ..
      //- double             DLARAN;
      //- Complex         ZLARND;
      // EXTERNAL DLARAN, ZLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, EXP, LOG, MAX, MIN, MOD

      if (N <= 0) return;
      zlaset('Full', N, N, CZERO, CZERO, A, LDA );

      // Insure a correct ISEED

      if( (ISEED( 4 ) % 2) != 1 ) ISEED( 4 ) = ISEED( 4 ) + 1;

      // Compute diagonal and subdiagonal according to ITYPE, NZ1, NZ2,
      // and RCOND

      if ( ITYPE != 0 ) {
         if ( ( ITYPE ).abs() >= 4 ) {
            KBEG = max( 1, min( N, NZ1+1 ) );
            KEND = max( KBEG, min( N, N-NZ2 ) );
            KLEN = KEND + 1 - KBEG;
         } else {
            KBEG = 1;
            KEND = N;
            KLEN = N;
         }
         ISDB = 1;
         ISDE = 0;
         GO TO ( 10, 30, 50, 80, 100, 120, 140, 160, 180, 200 )( ITYPE ).abs();

         // abs(ITYPE) = 1: Identity

         } // 10
         for (JD = 1; JD <= N; JD++) { // 20
            A[JD][JD] = CONE;
         } // 20
         GO TO 220;

         // abs(ITYPE) = 2: Transposed Jordan block

         } // 30
         for (JD = 1; JD <= N - 1; JD++) { // 40
            A[JD+1][JD] = CONE;
         } // 40
         ISDB = 1;
         ISDE = N - 1;
         GO TO 220;

         // abs(ITYPE) = 3: Transposed Jordan block, followed by the
                         // identity.

         } // 50
         K = ( N-1 ) / 2;
         for (JD = 1; JD <= K; JD++) { // 60
            A[JD+1][JD] = CONE;
         } // 60
         ISDB = 1;
         ISDE = K;
         for (JD = K + 2; JD <= 2*K + 1; JD++) { // 70
            A[JD][JD] = CONE;
         } // 70
         GO TO 220;

         // abs(ITYPE) = 4: 1,...,k

         } // 80
         for (JD = KBEG; JD <= KEND; JD++) { // 90
            A[JD][JD] = DCMPLX( JD-NZ1 );
         } // 90
         GO TO 220;

         // abs(ITYPE) = 5: One large D value:

         } // 100
         for (JD = KBEG + 1; JD <= KEND; JD++) { // 110
            A[JD][JD] = DCMPLX( RCOND );
         } // 110
         A[KBEG][KBEG] = CONE;
         GO TO 220;

         // abs(ITYPE) = 6: One small D value:

         } // 120
         for (JD = KBEG; JD <= KEND - 1; JD++) { // 130
            A[JD][JD] = CONE;
         } // 130
         A[KEND][KEND] = DCMPLX( RCOND );
         GO TO 220;

         // abs(ITYPE) = 7: Exponentially distributed D values:

         } // 140
         A[KBEG][KBEG] = CONE;
         if ( KLEN > 1 ) {
            ALPHA = RCOND**( ONE / (KLEN-1).toDouble() );
            for (I = 2; I <= KLEN; I++) { // 150
               A[NZ1+I][NZ1+I] = DCMPLX( ALPHA**(I-1).toDouble() );
            } // 150
         }
         GO TO 220;

         // abs(ITYPE) = 8: Arithmetically distributed D values:

         } // 160
         A[KBEG][KBEG] = CONE;
         if ( KLEN > 1 ) {
            ALPHA = ( ONE-RCOND ) / (KLEN-1).toDouble();
            for (I = 2; I <= KLEN; I++) { // 170
               A[NZ1+I][NZ1+I] = DCMPLX( (KLEN-I).toDouble()*ALPHA+RCOND );
            } // 170
         }
         GO TO 220;

         // abs(ITYPE) = 9: Randomly distributed D values on ( RCOND, 1):

         } // 180
         ALPHA = LOG( RCOND );
         for (JD = KBEG; JD <= KEND; JD++) { // 190
            A[JD][JD] = EXP( ALPHA*dlaran( ISEED ) );
         } // 190
         GO TO 220;

         // abs(ITYPE) = 10: Randomly distributed D values from DIST

         } // 200
         for (JD = KBEG; JD <= KEND; JD++) { // 210
            A[JD][JD] = ZLARND( IDIST, ISEED );
         } // 210

         } // 220

         // Scale by AMAGN

         for (JD = KBEG; JD <= KEND; JD++) { // 230
            A[JD][JD] = AMAGN*(A( JD, JD )).toDouble();
         } // 230
         for (JD = ISDB; JD <= ISDE; JD++) { // 240
            A[JD+1][JD] = AMAGN*(A( JD+1, JD )).toDouble();
         } // 240

         // If RSIGN = true , assign random signs to diagonal and
         // subdiagonal

         if ( RSIGN ) {
            for (JD = KBEG; JD <= KEND; JD++) { // 250
               if ( (A( JD, JD )).toDouble() != ZERO ) {
                  CTEMP = ZLARND( 3, ISEED );
                  CTEMP = CTEMP / ( CTEMP ).abs();
                  A[JD][JD] = CTEMP*(A( JD, JD )).toDouble();
               }
            } // 250
            for (JD = ISDB; JD <= ISDE; JD++) { // 260
               if ( (A( JD+1, JD )).toDouble() != ZERO ) {
                  CTEMP = ZLARND( 3, ISEED );
                  CTEMP = CTEMP / ( CTEMP ).abs();
                  A[JD+1][JD] = CTEMP*(A( JD+1, JD )).toDouble();
               }
            } // 260
         }

         // Reverse if ITYPE < 0

         if ( ITYPE < 0 ) {
            for (JD = KBEG; JD <= ( KBEG+KEND-1 ) / 2; JD++) { // 270
               CTEMP = A( JD, JD );
               A[JD][JD] = A( KBEG+KEND-JD, KBEG+KEND-JD );
               A[KBEG+KEND-JD][KBEG+KEND-JD] = CTEMP;
            } // 270
            for (JD = 1; JD <= ( N-1 ) / 2; JD++) { // 280
               CTEMP = A( JD+1, JD );
               A[JD+1][JD] = A( N+1-JD, N-JD );
               A[N+1-JD][N-JD] = CTEMP;
            } // 280
         }

      }

      // Fill in upper triangle

      if ( TRIANG != ZERO ) {
         for (JC = 2; JC <= N; JC++) { // 300
            for (JR = 1; JR <= JC - 1; JR++) { // 290
               A[JR][JC] = TRIANG*ZLARND( IDIST, ISEED );
            } // 290
         } // 300
      }

      }
