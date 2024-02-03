      SUBROUTINE DSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, WAY
      int                INFO, LDA, N
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Local Scalars ..
      LOGICAL            UPPER, CONVERT
      int                I, IP, J
      DOUBLE PRECISION   TEMP
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5

      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYCONV', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      IF( UPPER ) THEN
*
*      A is UPPER
*
*      Convert A (A is upper)
*
*        Convert VALUE
*
         IF ( CONVERT ) THEN
            I=N
            E(1)=ZERO
            DO WHILE ( I .GT. 1 )
               IF( IPIV(I) .LT. 0 ) THEN
                  E(I)=A(I-1,I)
                  E(I-1)=ZERO
                  A(I-1,I)=ZERO
                  I=I-1
               ELSE
                  E(I)=ZERO
               ENDIF
               I=I-1
            END DO
*
*        Convert PERMUTATIONS
*
         I=N
         DO WHILE ( I .GE. 1 )
            IF( IPIV(I) .GT. 0) THEN
               IP=IPIV(I)
               IF( I .LT. N) THEN
                  DO 12 J= I+1,N
                    TEMP=A(IP,J)
                    A(IP,J)=A(I,J)
                    A(I,J)=TEMP
 12            CONTINUE
               ENDIF
            ELSE
              IP=-IPIV(I)
               IF( I .LT. N) THEN
             DO 13 J= I+1,N
                 TEMP=A(IP,J)
                 A(IP,J)=A(I-1,J)
                 A(I-1,J)=TEMP
 13            CONTINUE
                ENDIF
                I=I-1
           ENDIF
           I=I-1
        END DO

         ELSE
*
*      Revert A (A is upper)
*
*
*        Revert PERMUTATIONS
*
            I=1
            DO WHILE ( I .LE. N )
               IF( IPIV(I) .GT. 0 ) THEN
                  IP=IPIV(I)
                  IF( I .LT. N) THEN
                  DO J= I+1,N
                    TEMP=A(IP,J)
                    A(IP,J)=A(I,J)
                    A(I,J)=TEMP
                  END DO
                  ENDIF
               ELSE
                 IP=-IPIV(I)
                 I=I+1
                 IF( I .LT. N) THEN
                    DO J= I+1,N
                       TEMP=A(IP,J)
                       A(IP,J)=A(I-1,J)
                       A(I-1,J)=TEMP
                    END DO
                 ENDIF
               ENDIF
               I=I+1
            END DO
*
*        Revert VALUE
*
            I=N
            DO WHILE ( I .GT. 1 )
               IF( IPIV(I) .LT. 0 ) THEN
                  A(I-1,I)=E(I)
                  I=I-1
               ENDIF
               I=I-1
            END DO
         END IF
      ELSE
*
*      A is LOWER
*
         IF ( CONVERT ) THEN
*
*      Convert A (A is lower)
*
*
*        Convert VALUE
*
            I=1
            E(N)=ZERO
            DO WHILE ( I .LE. N )
               IF( I.LT.N .AND. IPIV(I) .LT. 0 ) THEN
                  E(I)=A(I+1,I)
                  E(I+1)=ZERO
                  A(I+1,I)=ZERO
                  I=I+1
               ELSE
                  E(I)=ZERO
               ENDIF
               I=I+1
            END DO
*
*        Convert PERMUTATIONS
*
         I=1
         DO WHILE ( I .LE. N )
            IF( IPIV(I) .GT. 0 ) THEN
               IP=IPIV(I)
               IF (I .GT. 1) THEN
               DO 22 J= 1,I-1
                 TEMP=A(IP,J)
                 A(IP,J)=A(I,J)
                 A(I,J)=TEMP
 22            CONTINUE
               ENDIF
            ELSE
              IP=-IPIV(I)
              IF (I .GT. 1) THEN
              DO 23 J= 1,I-1
                 TEMP=A(IP,J)
                 A(IP,J)=A(I+1,J)
                 A(I+1,J)=TEMP
 23           CONTINUE
              ENDIF
              I=I+1
           ENDIF
           I=I+1
        END DO
         ELSE
*
*      Revert A (A is lower)
*
*
*        Revert PERMUTATIONS
*
            I=N
            DO WHILE ( I .GE. 1 )
               IF( IPIV(I) .GT. 0 ) THEN
                  IP=IPIV(I)
                  IF (I .GT. 1) THEN
                     DO J= 1,I-1
                        TEMP=A(I,J)
                        A(I,J)=A(IP,J)
                        A(IP,J)=TEMP
                     END DO
                  ENDIF
               ELSE
                  IP=-IPIV(I)
                  I=I-1
                  IF (I .GT. 1) THEN
                     DO J= 1,I-1
                        TEMP=A(I+1,J)
                        A(I+1,J)=A(IP,J)
                        A(IP,J)=TEMP
                     END DO
                  ENDIF
               ENDIF
               I=I-1
            END DO
*
*        Revert VALUE
*
            I=1
            DO WHILE ( I .LE. N-1 )
               IF( IPIV(I) .LT. 0 ) THEN
                  A(I+1,I)=E(I)
                  I=I+1
               ENDIF
               I=I+1
            END DO
         END IF
      END IF

      RETURN
*
*     End of DSYCONV
*
      END
