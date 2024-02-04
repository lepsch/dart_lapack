      void dlatm4(ITYPE, N, NZ1, NZ2, ISIGN, AMAGN, RCOND, TRIANG, IDIST, ISEED, A, LDA ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IDIST, ISIGN, ITYPE, LDA, N, NZ1, NZ2;
      double             AMAGN, RCOND, TRIANG;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      double             HALF;
      const              HALF = ONE / TWO ;
      // ..
      // .. Local Scalars ..
      int                I, IOFF, ISDB, ISDE, JC, JD, JR, K, KBEG, KEND, KLEN;
      double             ALPHA, CL, CR, SAFMIN, SL, SR, SV1, SV2, TEMP;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DLARAN, DLARND;
      // EXTERNAL DLAMCH, DLARAN, DLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, EXP, LOG, MAX, MIN, MOD, SQRT
      // ..
      // .. Executable Statements ..

      if (N <= 0) return;
      dlaset('Full', N, N, ZERO, ZERO, A, LDA );

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
            A[JD, JD] = ONE;
         } // 20
         GO TO 220;

         // abs(ITYPE) = 2: Transposed Jordan block

         } // 30
         for (JD = 1; JD <= N - 1; JD++) { // 40
            A[JD+1, JD] = ONE;
         } // 40
         ISDB = 1;
         ISDE = N - 1;
         GO TO 220;

         // abs(ITYPE) = 3: Transposed Jordan block, followed by the
                         // identity.

         } // 50
         K = ( N-1 ) / 2;
         for (JD = 1; JD <= K; JD++) { // 60
            A[JD+1, JD] = ONE;
         } // 60
         ISDB = 1;
         ISDE = K;
         for (JD = K + 2; JD <= 2*K + 1; JD++) { // 70
            A[JD, JD] = ONE;
         } // 70
         GO TO 220;

         // abs(ITYPE) = 4: 1,...,k

         } // 80
         for (JD = KBEG; JD <= KEND; JD++) { // 90
            A[JD, JD] = DBLE( JD-NZ1 );
         } // 90
         GO TO 220;

         // abs(ITYPE) = 5: One large D value:

         } // 100
         for (JD = KBEG + 1; JD <= KEND; JD++) { // 110
            A[JD, JD] = RCOND;
         } // 110
         A[KBEG, KBEG] = ONE;
         GO TO 220;

         // abs(ITYPE) = 6: One small D value:

         } // 120
         for (JD = KBEG; JD <= KEND - 1; JD++) { // 130
            A[JD, JD] = ONE;
         } // 130
         A[KEND, KEND] = RCOND;
         GO TO 220;

         // abs(ITYPE) = 7: Exponentially distributed D values:

         } // 140
         A[KBEG, KBEG] = ONE;
         if ( KLEN > 1 ) {
            ALPHA = RCOND**( ONE / DBLE( KLEN-1 ) );
            for (I = 2; I <= KLEN; I++) { // 150
               A[NZ1+I, NZ1+I] = ALPHA**DBLE( I-1 );
            } // 150
         }
         GO TO 220;

         // abs(ITYPE) = 8: Arithmetically distributed D values:

         } // 160
         A[KBEG, KBEG] = ONE;
         if ( KLEN > 1 ) {
            ALPHA = ( ONE-RCOND ) / DBLE( KLEN-1 );
            for (I = 2; I <= KLEN; I++) { // 170
               A[NZ1+I, NZ1+I] = DBLE( KLEN-I )*ALPHA + RCOND;
            } // 170
         }
         GO TO 220;

         // abs(ITYPE) = 9: Randomly distributed D values on ( RCOND, 1):

         } // 180
         ALPHA = LOG( RCOND );
         for (JD = KBEG; JD <= KEND; JD++) { // 190
            A[JD, JD] = EXP( ALPHA*DLARAN( ISEED ) );
         } // 190
         GO TO 220;

         // abs(ITYPE) = 10: Randomly distributed D values from DIST

         } // 200
         for (JD = KBEG; JD <= KEND; JD++) { // 210
            A[JD, JD] = DLARND( IDIST, ISEED );
         } // 210

         } // 220

         // Scale by AMAGN

         for (JD = KBEG; JD <= KEND; JD++) { // 230
            A[JD, JD] = AMAGN*DBLE( A( JD, JD ) );
         } // 230
         for (JD = ISDB; JD <= ISDE; JD++) { // 240
            A[JD+1, JD] = AMAGN*DBLE( A( JD+1, JD ) );
         } // 240

         // If ISIGN = 1 or 2, assign random signs to diagonal and
         // subdiagonal

         if ( ISIGN > 0 ) {
            for (JD = KBEG; JD <= KEND; JD++) { // 250
               if ( DBLE( A( JD, JD ) ) != ZERO ) {
                  if[DLARAN( ISEED ) > HALF ) A( JD, JD] = -A( JD, JD );
               }
            } // 250
            for (JD = ISDB; JD <= ISDE; JD++) { // 260
               if ( DBLE( A( JD+1, JD ) ) != ZERO ) {
                  if[DLARAN( ISEED ) > HALF ) A( JD+1, JD] = -A( JD+1, JD );
               }
            } // 260
         }

         // Reverse if ITYPE < 0

         if ( ITYPE < 0 ) {
            for (JD = KBEG; JD <= ( KBEG+KEND-1 ) / 2; JD++) { // 270
               TEMP = A( JD, JD );
               A[JD, JD] = A( KBEG+KEND-JD, KBEG+KEND-JD );
               A[KBEG+KEND-JD, KBEG+KEND-JD] = TEMP;
            } // 270
            for (JD = 1; JD <= ( N-1 ) / 2; JD++) { // 280
               TEMP = A( JD+1, JD );
               A[JD+1, JD] = A( N+1-JD, N-JD );
               A[N+1-JD, N-JD] = TEMP;
            } // 280
         }

         // If ISIGN = 2, and no subdiagonals already, then apply
         // random rotations to make 2x2 blocks.

         if ( ISIGN == 2 && ITYPE != 2 && ITYPE != 3 ) {
            SAFMIN = DLAMCH( 'S' );
            for (JD = KBEG; 2 < 0 ? JD >= KEND - 1 : JD <= KEND - 1; JD += 2) { // 290
               if ( DLARAN( ISEED ) > HALF ) {

                  // Rotation on left.

                  CL = TWO*DLARAN( ISEED ) - ONE;
                  SL = TWO*DLARAN( ISEED ) - ONE;
                  TEMP = ONE / max( SAFMIN, sqrt( CL**2+SL**2 ) );
                  CL = CL*TEMP;
                  SL = SL*TEMP;

                  // Rotation on right.

                  CR = TWO*DLARAN( ISEED ) - ONE;
                  SR = TWO*DLARAN( ISEED ) - ONE;
                  TEMP = ONE / max( SAFMIN, sqrt( CR**2+SR**2 ) );
                  CR = CR*TEMP;
                  SR = SR*TEMP;

                  // Apply

                  SV1 = A( JD, JD );
                  SV2 = A( JD+1, JD+1 );
                  A[JD, JD] = CL*CR*SV1 + SL*SR*SV2;
                  A[JD+1, JD] = -SL*CR*SV1 + CL*SR*SV2;
                  A[JD, JD+1] = -CL*SR*SV1 + SL*CR*SV2;
                  A[JD+1, JD+1] = SL*SR*SV1 + CL*CR*SV2;
               }
            } // 290
         }

      }

      // Fill in upper triangle (except for 2x2 blocks)

      if ( TRIANG != ZERO ) {
         if ( ISIGN != 2 || ITYPE == 2 || ITYPE == 3 ) {
            IOFF = 1;
         } else {
            IOFF = 2;
            for (JR = 1; JR <= N - 1; JR++) { // 300
               if( A( JR+1, JR ) == ZERO ) A( JR, JR+1 ) = TRIANG*DLARND( IDIST, ISEED );
            } // 300
         }

         for (JC = 2; JC <= N; JC++) { // 320
            for (JR = 1; JR <= JC - IOFF; JR++) { // 310
               A[JR, JC] = TRIANG*DLARND( IDIST, ISEED );
            } // 310
         } // 320
      }

      return;
      }
