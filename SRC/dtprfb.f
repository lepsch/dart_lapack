      SUBROUTINE DTPRFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V, LDV, T, LDT, A, LDA, B, LDB, WORK, LDWORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String    DIRECT, SIDE, STOREV, TRANS;
      int       K, L, LDA, LDB, LDT, LDV, LDWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), T( LDT, * ), V( LDV, * ), WORK( LDWORK, * );
      // ..

*  ==========================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const     ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int       I, J, MP, NP, KP;
      bool      LEFT, FORWARD, COLUMN, RIGHT, BACKWARD, ROW;
      // ..
      // .. External Functions ..
      bool      LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DTRMM
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( M.LE.0 .OR. N.LE.0 .OR. K.LE.0 .OR. L.LT.0 ) RETURN

      if ( LSAME( STOREV, 'C' ) ) {
         COLUMN = .TRUE.
         ROW = .FALSE.
      } else if ( LSAME( STOREV, 'R' ) ) {
         COLUMN = .FALSE.
         ROW = .TRUE.
      } else {
         COLUMN = .FALSE.
         ROW = .FALSE.
      }

      if ( LSAME( SIDE, 'L' ) ) {
         LEFT = .TRUE.
         RIGHT = .FALSE.
      } else if ( LSAME( SIDE, 'R' ) ) {
         LEFT = .FALSE.
         RIGHT = .TRUE.
      } else {
         LEFT = .FALSE.
         RIGHT = .FALSE.
      }

      if ( LSAME( DIRECT, 'F' ) ) {
         FORWARD = .TRUE.
         BACKWARD = .FALSE.
      } else if ( LSAME( DIRECT, 'B' ) ) {
         FORWARD = .FALSE.
         BACKWARD = .TRUE.
      } else {
         FORWARD = .FALSE.
         BACKWARD = .FALSE.
      }

* ---------------------------------------------------------------------------

      if ( COLUMN .AND. FORWARD .AND. LEFT  ) {

* ---------------------------------------------------------------------------

         // Let  W =  [ I ]    (K-by-K)
                   // [ V ]    (M-by-K)

         // Form  H C  or  H**T C  where  C = [ A ]  (K-by-N)
                                           // [ B ]  (M-by-N)

         // H = I - W T W**T          or  H**T = I - W T**T W**T

         // A = A -   T (A + V**T B)  or  A = A -   T**T (A + V**T B)
         // B = B - V T (A + V**T B)  or  B = B - V T**T (A + V**T B)

* ---------------------------------------------------------------------------

         MP = MIN( M-L+1, M )
         KP = MIN( L+1, K )

         DO J = 1, N
            DO I = 1, L
               WORK( I, J ) = B( M-L+I, J )
            END DO
         END DO
         dtrmm('L', 'U', 'T', 'N', L, N, ONE, V( MP, 1 ), LDV, WORK, LDWORK )          CALL DGEMM( 'T', 'N', L, N, M-L, ONE, V, LDV, B, LDB, ONE, WORK, LDWORK )          CALL DGEMM( 'T', 'N', K-L, N, M, ONE, V( 1, KP ), LDV, B, LDB, ZERO, WORK( KP, 1 ), LDWORK );

         DO J = 1, N
            DO I = 1, K
               WORK( I, J ) = WORK( I, J ) + A( I, J )
            END DO
         END DO

         dtrmm('L', 'U', TRANS, 'N', K, N, ONE, T, LDT, WORK, LDWORK );

         DO J = 1, N
            DO I = 1, K
               A( I, J ) = A( I, J ) - WORK( I, J )
            END DO
         END DO

         dgemm('N', 'N', M-L, N, K, -ONE, V, LDV, WORK, LDWORK, ONE, B, LDB )          CALL DGEMM( 'N', 'N', L, N, K-L, -ONE, V( MP, KP ), LDV, WORK( KP, 1 ), LDWORK, ONE, B( MP, 1 ),  LDB )          CALL DTRMM( 'L', 'U', 'N', 'N', L, N, ONE, V( MP, 1 ), LDV, WORK, LDWORK );
         DO J = 1, N
            DO I = 1, L
               B( M-L+I, J ) = B( M-L+I, J ) - WORK( I, J )
            END DO
         END DO

* ---------------------------------------------------------------------------

      } else if ( COLUMN .AND. FORWARD .AND. RIGHT ) {

* ---------------------------------------------------------------------------

         // Let  W =  [ I ]    (K-by-K)
                   // [ V ]    (N-by-K)

         // Form  C H or  C H**T  where  C = [ A B ] (A is M-by-K, B is M-by-N)

         // H = I - W T W**T          or  H**T = I - W T**T W**T

         // A = A - (A + B V) T      or  A = A - (A + B V) T**T
         // B = B - (A + B V) T V**T  or  B = B - (A + B V) T**T V**T

* ---------------------------------------------------------------------------

         NP = MIN( N-L+1, N )
         KP = MIN( L+1, K )

         DO J = 1, L
            DO I = 1, M
               WORK( I, J ) = B( I, N-L+J )
            END DO
         END DO
         dtrmm('R', 'U', 'N', 'N', M, L, ONE, V( NP, 1 ), LDV, WORK, LDWORK )          CALL DGEMM( 'N', 'N', M, L, N-L, ONE, B, LDB, V, LDV, ONE, WORK, LDWORK )          CALL DGEMM( 'N', 'N', M, K-L, N, ONE, B, LDB, V( 1, KP ), LDV, ZERO, WORK( 1, KP ), LDWORK );

         DO J = 1, K
            DO I = 1, M
               WORK( I, J ) = WORK( I, J ) + A( I, J )
            END DO
         END DO

         dtrmm('R', 'U', TRANS, 'N', M, K, ONE, T, LDT, WORK, LDWORK );

         DO J = 1, K
            DO I = 1, M
               A( I, J ) = A( I, J ) - WORK( I, J )
            END DO
         END DO

         dgemm('N', 'T', M, N-L, K, -ONE, WORK, LDWORK, V, LDV, ONE, B, LDB )          CALL DGEMM( 'N', 'T', M, L, K-L, -ONE, WORK( 1, KP ), LDWORK, V( NP, KP ), LDV, ONE, B( 1, NP ), LDB )          CALL DTRMM( 'R', 'U', 'T', 'N', M, L, ONE, V( NP, 1 ), LDV, WORK, LDWORK );
         DO J = 1, L
            DO I = 1, M
               B( I, N-L+J ) = B( I, N-L+J ) - WORK( I, J )
            END DO
         END DO

* ---------------------------------------------------------------------------

      } else if ( COLUMN .AND. BACKWARD .AND. LEFT ) {

* ---------------------------------------------------------------------------

         // Let  W =  [ V ]    (M-by-K)
                   // [ I ]    (K-by-K)

         // Form  H C  or  H**T C  where  C = [ B ]  (M-by-N)
                                           // [ A ]  (K-by-N)

         // H = I - W T W**T          or  H**T = I - W T**T W**T

         // A = A -   T (A + V**T B)  or  A = A -   T**T (A + V**T B)
         // B = B - V T (A + V**T B)  or  B = B - V T**T (A + V**T B)

* ---------------------------------------------------------------------------

         MP = MIN( L+1, M )
         KP = MIN( K-L+1, K )

         DO J = 1, N
            DO I = 1, L
               WORK( K-L+I, J ) = B( I, J )
            END DO
         END DO

         dtrmm('L', 'L', 'T', 'N', L, N, ONE, V( 1, KP ), LDV, WORK( KP, 1 ), LDWORK )          CALL DGEMM( 'T', 'N', L, N, M-L, ONE, V( MP, KP ), LDV, B( MP, 1 ), LDB, ONE, WORK( KP, 1 ), LDWORK )          CALL DGEMM( 'T', 'N', K-L, N, M, ONE, V, LDV, B, LDB, ZERO, WORK, LDWORK );

         DO J = 1, N
            DO I = 1, K
               WORK( I, J ) = WORK( I, J ) + A( I, J )
            END DO
         END DO

         dtrmm('L', 'L', TRANS, 'N', K, N, ONE, T, LDT, WORK, LDWORK );

         DO J = 1, N
            DO I = 1, K
               A( I, J ) = A( I, J ) - WORK( I, J )
            END DO
         END DO

         dgemm('N', 'N', M-L, N, K, -ONE, V( MP, 1 ), LDV, WORK, LDWORK, ONE, B( MP, 1 ), LDB )          CALL DGEMM( 'N', 'N', L, N, K-L, -ONE, V, LDV, WORK, LDWORK, ONE, B,  LDB )          CALL DTRMM( 'L', 'L', 'N', 'N', L, N, ONE, V( 1, KP ), LDV, WORK( KP, 1 ), LDWORK );
         DO J = 1, N
            DO I = 1, L
               B( I, J ) = B( I, J ) - WORK( K-L+I, J )
            END DO
         END DO

* ---------------------------------------------------------------------------

      } else if ( COLUMN .AND. BACKWARD .AND. RIGHT ) {

* ---------------------------------------------------------------------------

         // Let  W =  [ V ]    (N-by-K)
                   // [ I ]    (K-by-K)

         // Form  C H  or  C H**T  where  C = [ B A ] (B is M-by-N, A is M-by-K)

         // H = I - W T W**T          or  H**T = I - W T**T W**T

         // A = A - (A + B V) T      or  A = A - (A + B V) T**T
         // B = B - (A + B V) T V**T  or  B = B - (A + B V) T**T V**T

* ---------------------------------------------------------------------------

         NP = MIN( L+1, N )
         KP = MIN( K-L+1, K )

         DO J = 1, L
            DO I = 1, M
               WORK( I, K-L+J ) = B( I, J )
            END DO
         END DO
         dtrmm('R', 'L', 'N', 'N', M, L, ONE, V( 1, KP ), LDV, WORK( 1, KP ), LDWORK )          CALL DGEMM( 'N', 'N', M, L, N-L, ONE, B( 1, NP ), LDB, V( NP, KP ), LDV, ONE, WORK( 1, KP ), LDWORK )          CALL DGEMM( 'N', 'N', M, K-L, N, ONE, B, LDB, V, LDV, ZERO, WORK, LDWORK );

         DO J = 1, K
            DO I = 1, M
               WORK( I, J ) = WORK( I, J ) + A( I, J )
            END DO
         END DO

         dtrmm('R', 'L', TRANS, 'N', M, K, ONE, T, LDT, WORK, LDWORK );

         DO J = 1, K
            DO I = 1, M
               A( I, J ) = A( I, J ) - WORK( I, J )
            END DO
         END DO

         dgemm('N', 'T', M, N-L, K, -ONE, WORK, LDWORK, V( NP, 1 ), LDV, ONE, B( 1, NP ), LDB )          CALL DGEMM( 'N', 'T', M, L, K-L, -ONE, WORK, LDWORK, V, LDV, ONE, B, LDB )          CALL DTRMM( 'R', 'L', 'T', 'N', M, L, ONE, V( 1, KP ), LDV, WORK( 1, KP ), LDWORK );
         DO J = 1, L
            DO I = 1, M
               B( I, J ) = B( I, J ) - WORK( I, K-L+J )
            END DO
         END DO

* ---------------------------------------------------------------------------

      } else if ( ROW .AND. FORWARD .AND. LEFT ) {

* ---------------------------------------------------------------------------

         // Let  W =  [ I V ] ( I is K-by-K, V is K-by-M )

         // Form  H C  or  H**T C  where  C = [ A ]  (K-by-N)
                                           // [ B ]  (M-by-N)

         // H = I - W**T T W          or  H**T = I - W**T T**T W

         // A = A -     T (A + V B)  or  A = A -     T**T (A + V B)
         // B = B - V**T T (A + V B)  or  B = B - V**T T**T (A + V B)

* ---------------------------------------------------------------------------

         MP = MIN( M-L+1, M )
         KP = MIN( L+1, K )

         DO J = 1, N
            DO I = 1, L
               WORK( I, J ) = B( M-L+I, J )
            END DO
         END DO
         dtrmm('L', 'L', 'N', 'N', L, N, ONE, V( 1, MP ), LDV, WORK, LDB )          CALL DGEMM( 'N', 'N', L, N, M-L, ONE, V, LDV,B, LDB, ONE, WORK, LDWORK )          CALL DGEMM( 'N', 'N', K-L, N, M, ONE, V( KP, 1 ), LDV, B, LDB, ZERO, WORK( KP, 1 ), LDWORK );

         DO J = 1, N
            DO I = 1, K
               WORK( I, J ) = WORK( I, J ) + A( I, J )
            END DO
         END DO

         dtrmm('L', 'U', TRANS, 'N', K, N, ONE, T, LDT, WORK, LDWORK );

         DO J = 1, N
            DO I = 1, K
               A( I, J ) = A( I, J ) - WORK( I, J )
            END DO
         END DO

         dgemm('T', 'N', M-L, N, K, -ONE, V, LDV, WORK, LDWORK, ONE, B, LDB )          CALL DGEMM( 'T', 'N', L, N, K-L, -ONE, V( KP, MP ), LDV, WORK( KP, 1 ), LDWORK, ONE, B( MP, 1 ), LDB )          CALL DTRMM( 'L', 'L', 'T', 'N', L, N, ONE, V( 1, MP ), LDV, WORK, LDWORK );
         DO J = 1, N
            DO I = 1, L
               B( M-L+I, J ) = B( M-L+I, J ) - WORK( I, J )
            END DO
         END DO

* ---------------------------------------------------------------------------

      } else if ( ROW .AND. FORWARD .AND. RIGHT ) {

* ---------------------------------------------------------------------------

         // Let  W =  [ I V ] ( I is K-by-K, V is K-by-N )

         // Form  C H  or  C H**T  where  C = [ A B ] (A is M-by-K, B is M-by-N)

         // H = I - W**T T W            or  H**T = I - W**T T**T W

         // A = A - (A + B V**T) T      or  A = A - (A + B V**T) T**T
         // B = B - (A + B V**T) T V    or  B = B - (A + B V**T) T**T V

* ---------------------------------------------------------------------------

         NP = MIN( N-L+1, N )
         KP = MIN( L+1, K )

         DO J = 1, L
            DO I = 1, M
               WORK( I, J ) = B( I, N-L+J )
            END DO
         END DO
         dtrmm('R', 'L', 'T', 'N', M, L, ONE, V( 1, NP ), LDV, WORK, LDWORK )          CALL DGEMM( 'N', 'T', M, L, N-L, ONE, B, LDB, V, LDV, ONE, WORK, LDWORK )          CALL DGEMM( 'N', 'T', M, K-L, N, ONE, B, LDB, V( KP, 1 ), LDV, ZERO, WORK( 1, KP ), LDWORK );

         DO J = 1, K
            DO I = 1, M
               WORK( I, J ) = WORK( I, J ) + A( I, J )
            END DO
         END DO

         dtrmm('R', 'U', TRANS, 'N', M, K, ONE, T, LDT, WORK, LDWORK );

         DO J = 1, K
            DO I = 1, M
               A( I, J ) = A( I, J ) - WORK( I, J )
            END DO
         END DO

         dgemm('N', 'N', M, N-L, K, -ONE, WORK, LDWORK, V, LDV, ONE, B, LDB )          CALL DGEMM( 'N', 'N', M, L, K-L, -ONE, WORK( 1, KP ), LDWORK, V( KP, NP ), LDV, ONE, B( 1, NP ), LDB )          CALL DTRMM( 'R', 'L', 'N', 'N', M, L, ONE, V( 1, NP ), LDV, WORK, LDWORK );
         DO J = 1, L
            DO I = 1, M
               B( I, N-L+J ) = B( I, N-L+J ) - WORK( I, J )
            END DO
         END DO

* ---------------------------------------------------------------------------

      } else if ( ROW .AND. BACKWARD .AND. LEFT ) {

* ---------------------------------------------------------------------------

         // Let  W =  [ V I ] ( I is K-by-K, V is K-by-M )

         // Form  H C  or  H**T C  where  C = [ B ]  (M-by-N)
                                           // [ A ]  (K-by-N)

         // H = I - W**T T W          or  H**T = I - W**T T**T W

         // A = A -     T (A + V B)  or  A = A -     T**T (A + V B)
         // B = B - V**T T (A + V B)  or  B = B - V**T T**T (A + V B)

* ---------------------------------------------------------------------------

         MP = MIN( L+1, M )
         KP = MIN( K-L+1, K )

         DO J = 1, N
            DO I = 1, L
               WORK( K-L+I, J ) = B( I, J )
            END DO
         END DO
         dtrmm('L', 'U', 'N', 'N', L, N, ONE, V( KP, 1 ), LDV, WORK( KP, 1 ), LDWORK )          CALL DGEMM( 'N', 'N', L, N, M-L, ONE, V( KP, MP ), LDV, B( MP, 1 ), LDB, ONE, WORK( KP, 1 ), LDWORK )          CALL DGEMM( 'N', 'N', K-L, N, M, ONE, V, LDV, B, LDB, ZERO, WORK, LDWORK );

         DO J = 1, N
            DO I = 1, K
               WORK( I, J ) = WORK( I, J ) + A( I, J )
            END DO
         END DO

         dtrmm('L', 'L ', TRANS, 'N', K, N, ONE, T, LDT, WORK, LDWORK );

         DO J = 1, N
            DO I = 1, K
               A( I, J ) = A( I, J ) - WORK( I, J )
            END DO
         END DO

         dgemm('T', 'N', M-L, N, K, -ONE, V( 1, MP ), LDV, WORK, LDWORK, ONE, B( MP, 1 ), LDB )          CALL DGEMM( 'T', 'N', L, N, K-L, -ONE, V, LDV, WORK, LDWORK, ONE, B, LDB )          CALL DTRMM( 'L', 'U', 'T', 'N', L, N, ONE, V( KP, 1 ), LDV, WORK( KP, 1 ), LDWORK );
         DO J = 1, N
            DO I = 1, L
               B( I, J ) = B( I, J ) - WORK( K-L+I, J )
            END DO
         END DO

* ---------------------------------------------------------------------------

      } else if ( ROW .AND. BACKWARD .AND. RIGHT ) {

* ---------------------------------------------------------------------------

         // Let  W =  [ V I ] ( I is K-by-K, V is K-by-N )

         // Form  C H  or  C H**T  where  C = [ B A ] (A is M-by-K, B is M-by-N)

         // H = I - W**T T W            or  H**T = I - W**T T**T W

         // A = A - (A + B V**T) T      or  A = A - (A + B V**T) T**T
         // B = B - (A + B V**T) T V    or  B = B - (A + B V**T) T**T V

* ---------------------------------------------------------------------------

         NP = MIN( L+1, N )
         KP = MIN( K-L+1, K )

         DO J = 1, L
            DO I = 1, M
               WORK( I, K-L+J ) = B( I, J )
            END DO
         END DO
         dtrmm('R', 'U', 'T', 'N', M, L, ONE, V( KP, 1 ), LDV, WORK( 1, KP ), LDWORK )          CALL DGEMM( 'N', 'T', M, L, N-L, ONE, B( 1, NP ), LDB, V( KP, NP ), LDV, ONE, WORK( 1, KP ), LDWORK )          CALL DGEMM( 'N', 'T', M, K-L, N, ONE, B, LDB, V, LDV, ZERO, WORK, LDWORK );

         DO J = 1, K
            DO I = 1, M
               WORK( I, J ) = WORK( I, J ) + A( I, J )
            END DO
         END DO

         dtrmm('R', 'L', TRANS, 'N', M, K, ONE, T, LDT, WORK, LDWORK );

         DO J = 1, K
            DO I = 1, M
               A( I, J ) = A( I, J ) - WORK( I, J )
            END DO
         END DO

         dgemm('N', 'N', M, N-L, K, -ONE, WORK, LDWORK, V( 1, NP ), LDV, ONE, B( 1, NP ), LDB )          CALL DGEMM( 'N', 'N', M, L, K-L , -ONE, WORK, LDWORK, V, LDV, ONE, B, LDB )          CALL DTRMM( 'R', 'U', 'N', 'N', M, L, ONE, V( KP, 1 ), LDV, WORK( 1, KP ), LDWORK );
         DO J = 1, L
            DO I = 1, M
               B( I, J ) = B( I, J ) - WORK( I, K-L+J )
            END DO
         END DO

      }

      RETURN

      // End of DTPRFB

      }
