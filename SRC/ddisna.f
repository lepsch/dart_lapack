      SUBROUTINE DDISNA( JOB, M, N, D, SEP, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOB;
      int                INFO, M, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), SEP( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               DECR, EIGEN, INCR, LEFT, RIGHT, SING;
      int                I, K;
      double             ANORM, EPS, NEWGAP, OLDGAP, SAFMIN, THRESH;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH;
      // EXTERNAL LSAME, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      EIGEN = LSAME( JOB, 'E' )
      LEFT = LSAME( JOB, 'L' )
      RIGHT = LSAME( JOB, 'R' )
      SING = LEFT .OR. RIGHT
      if ( EIGEN ) {
         K = M
      } else if ( SING ) {
         K = MIN( M, N )
      }
      if ( .NOT.EIGEN && .NOT.SING ) {
         INFO = -1
      } else if ( M.LT.0 ) {
         INFO = -2
      } else if ( K.LT.0 ) {
         INFO = -3
      } else {
         INCR = true;
         DECR = true;
         for (I = 1; I <= K - 1; I++) { // 10
            if (INCR) INCR = INCR && D( I ).LE.D( I+1 )             IF( DECR ) DECR = DECR && D( I ).GE.D( I+1 );
         } // 10
         if ( SING && K.GT.0 ) {
            if (INCR) INCR = INCR && ZERO.LE.D( 1 )             IF( DECR ) DECR = DECR && D( K ).GE.ZERO;
         }
         IF( .NOT.( INCR .OR. DECR ) ) INFO = -4
      }
      if ( INFO != 0 ) {
         xerbla('DDISNA', -INFO );
         RETURN
      }

      // Quick return if possible

      if (K == 0) RETURN;

      // Compute reciprocal condition numbers

      if ( K == 1 ) {
         SEP( 1 ) = DLAMCH( 'O' )
      } else {
         OLDGAP = ABS( D( 2 )-D( 1 ) )
         SEP( 1 ) = OLDGAP
         for (I = 2; I <= K - 1; I++) { // 20
            NEWGAP = ABS( D( I+1 )-D( I ) )
            SEP( I ) = MIN( OLDGAP, NEWGAP )
            OLDGAP = NEWGAP
         } // 20
         SEP( K ) = OLDGAP
      }
      if ( SING ) {
         if ( ( LEFT && M.GT.N ) .OR. ( RIGHT && M.LT.N ) ) {
            if (INCR) SEP( 1 ) = MIN( SEP( 1 ), D( 1 ) )             IF( DECR ) SEP( K ) = MIN( SEP( K ), D( K ) );
         }
      }

      // Ensure that reciprocal condition numbers are not less than
      // threshold, in order to limit the size of the error bound

      EPS = DLAMCH( 'E' )
      SAFMIN = DLAMCH( 'S' )
      ANORM = MAX( ABS( D( 1 ) ), ABS( D( K ) ) )
      if ( ANORM == ZERO ) {
         THRESH = EPS
      } else {
         THRESH = MAX( EPS*ANORM, SAFMIN )
      }
      for (I = 1; I <= K; I++) { // 30
         SEP( I ) = MAX( SEP( I ), THRESH )
      } // 30

      RETURN

      // End of DDISNA

      }
