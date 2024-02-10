      void sdisna(final int JOB, final int M, final int N, final int D, final int SEP, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOB;
      int                INFO, M, N;
      double               D( * ), SEP( * );
      // ..

      double               ZERO;
      const              ZERO = 0.0 ;
      bool               DECR, EIGEN, INCR, LEFT, RIGHT, SING;
      int                I, K;
      double               ANORM, EPS, NEWGAP, OLDGAP, SAFMIN, THRESH;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH;
      // EXTERNAL lsame, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA

      // Test the input arguments

      INFO = 0;
      EIGEN = lsame( JOB, 'E' );
      LEFT = lsame( JOB, 'L' );
      RIGHT = lsame( JOB, 'R' );
      SING = LEFT || RIGHT;
      if ( EIGEN ) {
         K = M;
      } else if ( SING ) {
         K = min( M, N );
      }
      if ( !EIGEN && !SING ) {
         INFO = -1;
      } else if ( M < 0 ) {
         INFO = -2;
      } else if ( K < 0 ) {
         INFO = -3;
      } else {
         INCR = true;
         DECR = true;
         for (I = 1; I <= K - 1; I++) { // 10
            if (INCR) INCR = INCR && D( I ) <= D( I+1 );
            IF( DECR ) DECR = DECR && D( I ) >= D( I+1 );
         } // 10
         if ( SING && K > 0 ) {
            if (INCR) INCR = INCR && ZERO <= D( 1 );
            IF( DECR ) DECR = DECR && D( K ) >= ZERO;
         }
         if( !( INCR || DECR ) ) INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('SDISNA', -INFO );
         return;
      }

      // Quick return if possible

      if (K == 0) return;

      // Compute reciprocal condition numbers

      if ( K == 1 ) {
         SEP[1] = SLAMCH( 'O' );
      } else {
         OLDGAP = ABS( D( 2 )-D( 1 ) );
         SEP[1] = OLDGAP;
         for (I = 2; I <= K - 1; I++) { // 20
            NEWGAP = ABS( D( I+1 )-D( I ) );
            SEP[I] = min( OLDGAP, NEWGAP );
            OLDGAP = NEWGAP;
         } // 20
         SEP[K] = OLDGAP;
      }
      if ( SING ) {
         if ( ( LEFT && M > N ) || ( RIGHT && M < N ) ) {
            if (INCR) SEP( 1 ) = min( SEP( 1 ), D( 1 ) );
            IF[DECR ) SEP( K] = min( SEP( K ), D( K ) );
         }
      }

      // Ensure that reciprocal condition numbers are not less than
      // threshold, in order to limit the size of the error bound

      EPS = SLAMCH( 'E' );
      SAFMIN = SLAMCH( 'S' );
      ANORM = max( ( D( 1 ) ).abs(), ( D( K ) ).abs() );
      if ( ANORM == ZERO ) {
         THRESH = EPS;
      } else {
         THRESH = max( EPS*ANORM, SAFMIN );
      }
      for (I = 1; I <= K; I++) { // 30
         SEP[I] = max( SEP( I ), THRESH );
      } // 30

      }
