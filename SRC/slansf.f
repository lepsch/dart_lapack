      REAL FUNCTION SLANSF( NORM, TRANSR, UPLO, N, A, WORK )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          NORM, TRANSR, UPLO
      INTEGER            N
*     ..
*     .. Array Arguments ..
      REAL               A( 0: * ), WORK( 0: * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, IFM, ILU, NOE, N1, K, L, LDA
      REAL               SCALE, S, VALUE, AA, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, SISNAN
      EXTERNAL           LSAME, SISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLASSQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.EQ.0 ) THEN
         SLANSF = ZERO
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         SLANSF = ABS( A(0) )
         RETURN
      END IF
*
*     set noe = 1 if n is odd. if n is even set noe=0
*
      NOE = 1
      IF( MOD( N, 2 ).EQ.0 )
     $   NOE = 0
*
*     set ifm = 0 when form='T or 't' and 1 otherwise
*
      IFM = 1
      IF( LSAME( TRANSR, 'T' ) )
     $   IFM = 0
*
*     set ilu = 0 when uplo='U or 'u' and 1 otherwise
*
      ILU = 1
      IF( LSAME( UPLO, 'U' ) )
     $   ILU = 0
*
*     set lda = (n+1)/2 when ifm = 0
*     set lda = n when ifm = 1 and noe = 1
*     set lda = n+1 when ifm = 1 and noe = 0
*
      IF( IFM.EQ.1 ) THEN
         IF( NOE.EQ.1 ) THEN
            LDA = N
         ELSE
*           noe=0
            LDA = N + 1
         END IF
      ELSE
*        ifm=0
         LDA = ( N+1 ) / 2
      END IF
*
      IF( LSAME( NORM, 'M' ) ) THEN
*
*       Find max(abs(A(i,j))).
*
         K = ( N+1 ) / 2
         VALUE = ZERO
         IF( NOE.EQ.1 ) THEN
*           n is odd
            IF( IFM.EQ.1 ) THEN
*           A is n by k
               DO J = 0, K - 1
                  DO I = 0, N - 1
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) )
     $                    VALUE = TEMP
                  END DO
               END DO
            ELSE
*              xpose case; A is k by n
               DO J = 0, N - 1
                  DO I = 0, K - 1
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) )
     $                    VALUE = TEMP
                  END DO
               END DO
            END IF
         ELSE
*           n is even
            IF( IFM.EQ.1 ) THEN
*              A is n+1 by k
               DO J = 0, K - 1
                  DO I = 0, N
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) )
     $                    VALUE = TEMP
                  END DO
               END DO
            ELSE
*              xpose case; A is k by n+1
               DO J = 0, N
                  DO I = 0, K - 1
                     TEMP = ABS( A( I+J*LDA ) )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) )
     $                    VALUE = TEMP
                  END DO
               END DO
            END IF
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR.
     $         ( NORM.EQ.'1' ) ) THEN
*
*        Find normI(A) ( = norm1(A), since A is symmetric).
*
         IF( IFM.EQ.1 ) THEN
            K = N / 2
            IF( NOE.EQ.1 ) THEN
*              n is odd
               IF( ILU.EQ.0 ) THEN
                  DO I = 0, K - 1
                     WORK( I ) = ZERO
                  END DO
                  DO J = 0, K
                     S = ZERO
                     DO I = 0, K + J - 1
                        AA = ABS( A( I+J*LDA ) )
*                       -> A(i,j+k)
                        S = S + AA
                        WORK( I ) = WORK( I ) + AA
                     END DO
                     AA = ABS( A( I+J*LDA ) )
*                    -> A(j+k,j+k)
                     WORK( J+K ) = S + AA
                     IF( I.EQ.K+K )
     $                  GO TO 10
                     I = I + 1
                     AA = ABS( A( I+J*LDA ) )
*                    -> A(j,j)
                     WORK( J ) = WORK( J ) + AA
                     S = ZERO
                     DO L = J + 1, K - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
*                       -> A(l,j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     END DO
                     WORK( J ) = WORK( J ) + S
                  END DO
   10             CONTINUE
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) )
     $                    VALUE = TEMP
                  END DO
               ELSE
*                 ilu = 1
                  K = K + 1
*                 k=(n+1)/2 for n odd and ilu=1
                  DO I = K, N - 1
                     WORK( I ) = ZERO
                  END DO
                  DO J = K - 1, 0, -1
                     S = ZERO
                     DO I = 0, J - 2
                        AA = ABS( A( I+J*LDA ) )
*                       -> A(j+k,i+k)
                        S = S + AA
                        WORK( I+K ) = WORK( I+K ) + AA
                     END DO
                     IF( J.GT.0 ) THEN
                        AA = ABS( A( I+J*LDA ) )
*                       -> A(j+k,j+k)
                        S = S + AA
                        WORK( I+K ) = WORK( I+K ) + S
*                       i=j
                        I = I + 1
                     END IF
                     AA = ABS( A( I+J*LDA ) )
*                    -> A(j,j)
                     WORK( J ) = AA
                     S = ZERO
                     DO L = J + 1, N - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
*                       -> A(l,j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     END DO
                     WORK( J ) = WORK( J ) + S
                  END DO
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) )
     $                    VALUE = TEMP
                  END DO
               END IF
            ELSE
*              n is even
               IF( ILU.EQ.0 ) THEN
                  DO I = 0, K - 1
                     WORK( I ) = ZERO
                  END DO
                  DO J = 0, K - 1
                     S = ZERO
                     DO I = 0, K + J - 1
                        AA = ABS( A( I+J*LDA ) )
*                       -> A(i,j+k)
                        S = S + AA
                        WORK( I ) = WORK( I ) + AA
                     END DO
                     AA = ABS( A( I+J*LDA ) )
*                    -> A(j+k,j+k)
                     WORK( J+K ) = S + AA
                     I = I + 1
                     AA = ABS( A( I+J*LDA ) )
*                    -> A(j,j)
                     WORK( J ) = WORK( J ) + AA
                     S = ZERO
                     DO L = J + 1, K - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
*                       -> A(l,j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     END DO
                     WORK( J ) = WORK( J ) + S
                  END DO
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) )
     $                    VALUE = TEMP
                  END DO
               ELSE
*                 ilu = 1
                  DO I = K, N - 1
                     WORK( I ) = ZERO
                  END DO
                  DO J = K - 1, 0, -1
                     S = ZERO
                     DO I = 0, J - 1
                        AA = ABS( A( I+J*LDA ) )
*                       -> A(j+k,i+k)
                        S = S + AA
                        WORK( I+K ) = WORK( I+K ) + AA
                     END DO
                     AA = ABS( A( I+J*LDA ) )
*                    -> A(j+k,j+k)
                     S = S + AA
                     WORK( I+K ) = WORK( I+K ) + S
*                    i=j
                     I = I + 1
                     AA = ABS( A( I+J*LDA ) )
*                    -> A(j,j)
                     WORK( J ) = AA
                     S = ZERO
                     DO L = J + 1, N - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
*                       -> A(l,j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     END DO
                     WORK( J ) = WORK( J ) + S
                  END DO
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) )
     $                    VALUE = TEMP
                  END DO
               END IF
            END IF
         ELSE
*           ifm=0
            K = N / 2
            IF( NOE.EQ.1 ) THEN
*              n is odd
               IF( ILU.EQ.0 ) THEN
                  N1 = K
*                 n/2
                  K = K + 1
*                 k is the row size and lda
                  DO I = N1, N - 1
                     WORK( I ) = ZERO
                  END DO
                  DO J = 0, N1 - 1
                     S = ZERO
                     DO I = 0, K - 1
                        AA = ABS( A( I+J*LDA ) )
*                       A(j,n1+i)
                        WORK( I+N1 ) = WORK( I+N1 ) + AA
                        S = S + AA
                     END DO
                     WORK( J ) = S
                  END DO
*                 j=n1=k-1 is special
                  S = ABS( A( 0+J*LDA ) )
*                 A(k-1,k-1)
                  DO I = 1, K - 1
                     AA = ABS( A( I+J*LDA ) )
*                    A(k-1,i+n1)
                     WORK( I+N1 ) = WORK( I+N1 ) + AA
                     S = S + AA
                  END DO
                  WORK( J ) = WORK( J ) + S
                  DO J = K, N - 1
                     S = ZERO
                     DO I = 0, J - K - 1
                        AA = ABS( A( I+J*LDA ) )
*                       A(i,j-k)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     END DO
*                    i=j-k
                     AA = ABS( A( I+J*LDA ) )
*                    A(j-k,j-k)
                     S = S + AA
                     WORK( J-K ) = WORK( J-K ) + S
                     I = I + 1
                     S = ABS( A( I+J*LDA ) )
*                    A(j,j)
                     DO L = J + 1, N - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
*                       A(j,l)
                        WORK( L ) = WORK( L ) + AA
                        S = S + AA
                     END DO
                     WORK( J ) = WORK( J ) + S
                  END DO
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) )
     $                    VALUE = TEMP
                  END DO
               ELSE
*                 ilu=1
                  K = K + 1
*                 k=(n+1)/2 for n odd and ilu=1
                  DO I = K, N - 1
                     WORK( I ) = ZERO
                  END DO
                  DO J = 0, K - 2
*                    process
                     S = ZERO
                     DO I = 0, J - 1
                        AA = ABS( A( I+J*LDA ) )
*                       A(j,i)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     END DO
                     AA = ABS( A( I+J*LDA ) )
*                    i=j so process of A(j,j)
                     S = S + AA
                     WORK( J ) = S
*                    is initialised here
                     I = I + 1
*                    i=j process A(j+k,j+k)
                     AA = ABS( A( I+J*LDA ) )
                     S = AA
                     DO L = K + J + 1, N - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
*                       A(l,k+j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     END DO
                     WORK( K+J ) = WORK( K+J ) + S
                  END DO
*                 j=k-1 is special :process col A(k-1,0:k-1)
                  S = ZERO
                  DO I = 0, K - 2
                     AA = ABS( A( I+J*LDA ) )
*                    A(k,i)
                     WORK( I ) = WORK( I ) + AA
                     S = S + AA
                  END DO
*                 i=k-1
                  AA = ABS( A( I+J*LDA ) )
*                 A(k-1,k-1)
                  S = S + AA
                  WORK( I ) = S
*                 done with col j=k+1
                  DO J = K, N - 1
*                    process col j of A = A(j,0:k-1)
                     S = ZERO
                     DO I = 0, K - 1
                        AA = ABS( A( I+J*LDA ) )
*                       A(j,i)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     END DO
                     WORK( J ) = WORK( J ) + S
                  END DO
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) )
     $                    VALUE = TEMP
                  END DO
               END IF
            ELSE
*              n is even
               IF( ILU.EQ.0 ) THEN
                  DO I = K, N - 1
                     WORK( I ) = ZERO
                  END DO
                  DO J = 0, K - 1
                     S = ZERO
                     DO I = 0, K - 1
                        AA = ABS( A( I+J*LDA ) )
*                       A(j,i+k)
                        WORK( I+K ) = WORK( I+K ) + AA
                        S = S + AA
                     END DO
                     WORK( J ) = S
                  END DO
*                 j=k
                  AA = ABS( A( 0+J*LDA ) )
*                 A(k,k)
                  S = AA
                  DO I = 1, K - 1
                     AA = ABS( A( I+J*LDA ) )
*                    A(k,k+i)
                     WORK( I+K ) = WORK( I+K ) + AA
                     S = S + AA
                  END DO
                  WORK( J ) = WORK( J ) + S
                  DO J = K + 1, N - 1
                     S = ZERO
                     DO I = 0, J - 2 - K
                        AA = ABS( A( I+J*LDA ) )
*                       A(i,j-k-1)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     END DO
*                     i=j-1-k
                     AA = ABS( A( I+J*LDA ) )
*                    A(j-k-1,j-k-1)
                     S = S + AA
                     WORK( J-K-1 ) = WORK( J-K-1 ) + S
                     I = I + 1
                     AA = ABS( A( I+J*LDA ) )
*                    A(j,j)
                     S = AA
                     DO L = J + 1, N - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
*                       A(j,l)
                        WORK( L ) = WORK( L ) + AA
                        S = S + AA
                     END DO
                     WORK( J ) = WORK( J ) + S
                  END DO
*                 j=n
                  S = ZERO
                  DO I = 0, K - 2
                     AA = ABS( A( I+J*LDA ) )
*                    A(i,k-1)
                     WORK( I ) = WORK( I ) + AA
                     S = S + AA
                  END DO
*                 i=k-1
                  AA = ABS( A( I+J*LDA ) )
*                 A(k-1,k-1)
                  S = S + AA
                  WORK( I ) = WORK( I ) + S
                  VALUE = WORK ( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) )
     $                    VALUE = TEMP
                  END DO
               ELSE
*                 ilu=1
                  DO I = K, N - 1
                     WORK( I ) = ZERO
                  END DO
*                 j=0 is special :process col A(k:n-1,k)
                  S = ABS( A( 0 ) )
*                 A(k,k)
                  DO I = 1, K - 1
                     AA = ABS( A( I ) )
*                    A(k+i,k)
                     WORK( I+K ) = WORK( I+K ) + AA
                     S = S + AA
                  END DO
                  WORK( K ) = WORK( K ) + S
                  DO J = 1, K - 1
*                    process
                     S = ZERO
                     DO I = 0, J - 2
                        AA = ABS( A( I+J*LDA ) )
*                       A(j-1,i)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     END DO
                     AA = ABS( A( I+J*LDA ) )
*                    i=j-1 so process of A(j-1,j-1)
                     S = S + AA
                     WORK( J-1 ) = S
*                    is initialised here
                     I = I + 1
*                    i=j process A(j+k,j+k)
                     AA = ABS( A( I+J*LDA ) )
                     S = AA
                     DO L = K + J + 1, N - 1
                        I = I + 1
                        AA = ABS( A( I+J*LDA ) )
*                       A(l,k+j)
                        S = S + AA
                        WORK( L ) = WORK( L ) + AA
                     END DO
                     WORK( K+J ) = WORK( K+J ) + S
                  END DO
*                 j=k is special :process col A(k,0:k-1)
                  S = ZERO
                  DO I = 0, K - 2
                     AA = ABS( A( I+J*LDA ) )
*                    A(k,i)
                     WORK( I ) = WORK( I ) + AA
                     S = S + AA
                  END DO
*                 i=k-1
                  AA = ABS( A( I+J*LDA ) )
*                 A(k-1,k-1)
                  S = S + AA
                  WORK( I ) = S
*                 done with col j=k+1
                  DO J = K + 1, N
*                    process col j-1 of A = A(j-1,0:k-1)
                     S = ZERO
                     DO I = 0, K - 1
                        AA = ABS( A( I+J*LDA ) )
*                       A(j-1,i)
                        WORK( I ) = WORK( I ) + AA
                        S = S + AA
                     END DO
                     WORK( J-1 ) = WORK( J-1 ) + S
                  END DO
                  VALUE = WORK( 0 )
                  DO I = 1, N-1
                     TEMP = WORK( I )
                     IF( VALUE .LT. TEMP .OR. SISNAN( TEMP ) )
     $                    VALUE = TEMP
                  END DO
               END IF
            END IF
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*       Find normF(A).
*
         K = ( N+1 ) / 2
         SCALE = ZERO
         S = ONE
         IF( NOE.EQ.1 ) THEN
*           n is odd
            IF( IFM.EQ.1 ) THEN
*              A is normal
               IF( ILU.EQ.0 ) THEN
*                 A is upper
                  DO J = 0, K - 3
                     CALL SLASSQ( K-J-2, A( K+J+1+J*LDA ), 1, SCALE, S )
*                    L at A(k,0)
                  END DO
                  DO J = 0, K - 1
                     CALL SLASSQ( K+J-1, A( 0+J*LDA ), 1, SCALE, S )
*                    trap U at A(0,0)
                  END DO
                  S = S + S
*                 double s for the off diagonal elements
                  CALL SLASSQ( K-1, A( K ), LDA+1, SCALE, S )
*                 tri L at A(k,0)
                  CALL SLASSQ( K, A( K-1 ), LDA+1, SCALE, S )
*                 tri U at A(k-1,0)
               ELSE
*                 ilu=1 & A is lower
                  DO J = 0, K - 1
                     CALL SLASSQ( N-J-1, A( J+1+J*LDA ), 1, SCALE, S )
*                    trap L at A(0,0)
                  END DO
                  DO J = 0, K - 2
                     CALL SLASSQ( J, A( 0+( 1+J )*LDA ), 1, SCALE, S )
*                    U at A(0,1)
                  END DO
                  S = S + S
*                 double s for the off diagonal elements
                  CALL SLASSQ( K, A( 0 ), LDA+1, SCALE, S )
*                 tri L at A(0,0)
                  CALL SLASSQ( K-1, A( 0+LDA ), LDA+1, SCALE, S )
*                 tri U at A(0,1)
               END IF
            ELSE
*              A is xpose
               IF( ILU.EQ.0 ) THEN
*                 A**T is upper
                  DO J = 1, K - 2
                     CALL SLASSQ( J, A( 0+( K+J )*LDA ), 1, SCALE, S )
*                    U at A(0,k)
                  END DO
                  DO J = 0, K - 2
                     CALL SLASSQ( K, A( 0+J*LDA ), 1, SCALE, S )
*                    k by k-1 rect. at A(0,0)
                  END DO
                  DO J = 0, K - 2
                     CALL SLASSQ( K-J-1, A( J+1+( J+K-1 )*LDA ), 1,
     $                            SCALE, S )
*                    L at A(0,k-1)
                  END DO
                  S = S + S
*                 double s for the off diagonal elements
                  CALL SLASSQ( K-1, A( 0+K*LDA ), LDA+1, SCALE, S )
*                 tri U at A(0,k)
                  CALL SLASSQ( K, A( 0+( K-1 )*LDA ), LDA+1, SCALE, S )
*                 tri L at A(0,k-1)
               ELSE
*                 A**T is lower
                  DO J = 1, K - 1
                     CALL SLASSQ( J, A( 0+J*LDA ), 1, SCALE, S )
*                    U at A(0,0)
                  END DO
                  DO J = K, N - 1
                     CALL SLASSQ( K, A( 0+J*LDA ), 1, SCALE, S )
*                    k by k-1 rect. at A(0,k)
                  END DO
                  DO J = 0, K - 3
                     CALL SLASSQ( K-J-2, A( J+2+J*LDA ), 1, SCALE, S )
*                    L at A(1,0)
                  END DO
                  S = S + S
*                 double s for the off diagonal elements
                  CALL SLASSQ( K, A( 0 ), LDA+1, SCALE, S )
*                 tri U at A(0,0)
                  CALL SLASSQ( K-1, A( 1 ), LDA+1, SCALE, S )
*                 tri L at A(1,0)
               END IF
            END IF
         ELSE
*           n is even
            IF( IFM.EQ.1 ) THEN
*              A is normal
               IF( ILU.EQ.0 ) THEN
*                 A is upper
                  DO J = 0, K - 2
                     CALL SLASSQ( K-J-1, A( K+J+2+J*LDA ), 1, SCALE, S )
*                    L at A(k+1,0)
                  END DO
                  DO J = 0, K - 1
                     CALL SLASSQ( K+J, A( 0+J*LDA ), 1, SCALE, S )
*                    trap U at A(0,0)
                  END DO
                  S = S + S
*                 double s for the off diagonal elements
                  CALL SLASSQ( K, A( K+1 ), LDA+1, SCALE, S )
*                 tri L at A(k+1,0)
                  CALL SLASSQ( K, A( K ), LDA+1, SCALE, S )
*                 tri U at A(k,0)
               ELSE
*                 ilu=1 & A is lower
                  DO J = 0, K - 1
                     CALL SLASSQ( N-J-1, A( J+2+J*LDA ), 1, SCALE, S )
*                    trap L at A(1,0)
                  END DO
                  DO J = 1, K - 1
                     CALL SLASSQ( J, A( 0+J*LDA ), 1, SCALE, S )
*                    U at A(0,0)
                  END DO
                  S = S + S
*                 double s for the off diagonal elements
                  CALL SLASSQ( K, A( 1 ), LDA+1, SCALE, S )
*                 tri L at A(1,0)
                  CALL SLASSQ( K, A( 0 ), LDA+1, SCALE, S )
*                 tri U at A(0,0)
               END IF
            ELSE
*              A is xpose
               IF( ILU.EQ.0 ) THEN
*                 A**T is upper
                  DO J = 1, K - 1
                     CALL SLASSQ( J, A( 0+( K+1+J )*LDA ), 1, SCALE, S )
*                    U at A(0,k+1)
                  END DO
                  DO J = 0, K - 1
                     CALL SLASSQ( K, A( 0+J*LDA ), 1, SCALE, S )
*                    k by k rect. at A(0,0)
                  END DO
                  DO J = 0, K - 2
                     CALL SLASSQ( K-J-1, A( J+1+( J+K )*LDA ), 1, SCALE,
     $                            S )
*                    L at A(0,k)
                  END DO
                  S = S + S
*                 double s for the off diagonal elements
                  CALL SLASSQ( K, A( 0+( K+1 )*LDA ), LDA+1, SCALE, S )
*                 tri U at A(0,k+1)
                  CALL SLASSQ( K, A( 0+K*LDA ), LDA+1, SCALE, S )
*                 tri L at A(0,k)
               ELSE
*                 A**T is lower
                  DO J = 1, K - 1
                     CALL SLASSQ( J, A( 0+( J+1 )*LDA ), 1, SCALE, S )
*                    U at A(0,1)
                  END DO
                  DO J = K + 1, N
                     CALL SLASSQ( K, A( 0+J*LDA ), 1, SCALE, S )
*                    k by k rect. at A(0,k+1)
                  END DO
                  DO J = 0, K - 2
                     CALL SLASSQ( K-J-1, A( J+1+J*LDA ), 1, SCALE, S )
*                    L at A(0,0)
                  END DO
                  S = S + S
*                 double s for the off diagonal elements
                  CALL SLASSQ( K, A( LDA ), LDA+1, SCALE, S )
*                 tri L at A(0,1)
                  CALL SLASSQ( K, A( 0 ), LDA+1, SCALE, S )
*                 tri U at A(0,0)
               END IF
            END IF
         END IF
         VALUE = SCALE*SQRT( S )
      END IF
*
      SLANSF = VALUE
      RETURN
*
*     End of SLANSF
*
      END