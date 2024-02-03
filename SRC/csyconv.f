      SUBROUTINE CSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO, WAY;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), E( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO
      const              ZERO = (0.0E+0,0.0E+0) ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME

      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // .. Local Scalars ..
      bool               UPPER, CONVERT;
      int                I, IP, J;
      COMPLEX            TEMP
      // ..
      // .. Executable Statements ..

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5

      }
      if ( INFO.NE.0 ) {
         xerbla('CSYCONV', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      if ( UPPER ) {

       // A is UPPER

       // Convert A (A is upper)

         // Convert VALUE

         if ( CONVERT ) {
            I=N
            E(1)=ZERO
            DO WHILE ( I .GT. 1 )
               if ( IPIV(I) .LT. 0 ) {
                  E(I)=A(I-1,I)
                  E(I-1)=ZERO
                  A(I-1,I)=ZERO
                  I=I-1
               } else {
                  E(I)=ZERO
               ENDIF
               I=I-1
            END DO

         // Convert PERMUTATIONS

         I=N
         DO WHILE ( I .GE. 1 )
            if ( IPIV(I) .GT. 0) {
               IP=IPIV(I)
               if ( I .LT. N) {
                  DO 12 J= I+1,N
                    TEMP=A(IP,J)
                    A(IP,J)=A(I,J)
                    A(I,J)=TEMP
               } // 12
               ENDIF
            } else {
              IP=-IPIV(I)
               if ( I .LT. N) {
             DO 13 J= I+1,N
                 TEMP=A(IP,J)
                 A(IP,J)=A(I-1,J)
                 A(I-1,J)=TEMP
               } // 13
                ENDIF
                I=I-1
           ENDIF
           I=I-1
        END DO

         } else {

       // Revert A (A is upper)


         // Revert PERMUTATIONS

            I=1
            DO WHILE ( I .LE. N )
               if ( IPIV(I) .GT. 0 ) {
                  IP=IPIV(I)
                  if ( I .LT. N) {
                  DO J= I+1,N
                    TEMP=A(IP,J)
                    A(IP,J)=A(I,J)
                    A(I,J)=TEMP
                  END DO
                  ENDIF
               } else {
                 IP=-IPIV(I)
                 I=I+1
                 if ( I .LT. N) {
                    DO J= I+1,N
                       TEMP=A(IP,J)
                       A(IP,J)=A(I-1,J)
                       A(I-1,J)=TEMP
                    END DO
                 ENDIF
               ENDIF
               I=I+1
            END DO

         // Revert VALUE

            I=N
            DO WHILE ( I .GT. 1 )
               if ( IPIV(I) .LT. 0 ) {
                  A(I-1,I)=E(I)
                  I=I-1
               ENDIF
               I=I-1
            END DO
         }
      } else {

       // A is LOWER

         if ( CONVERT ) {

       // Convert A (A is lower)


         // Convert VALUE

            I=1
            E(N)=ZERO
            DO WHILE ( I .LE. N )
               if ( I.LT.N .AND. IPIV(I) .LT. 0 ) {
                  E(I)=A(I+1,I)
                  E(I+1)=ZERO
                  A(I+1,I)=ZERO
                  I=I+1
               } else {
                  E(I)=ZERO
               ENDIF
               I=I+1
            END DO

         // Convert PERMUTATIONS

         I=1
         DO WHILE ( I .LE. N )
            if ( IPIV(I) .GT. 0 ) {
               IP=IPIV(I)
               if (I .GT. 1) {
               DO 22 J= 1,I-1
                 TEMP=A(IP,J)
                 A(IP,J)=A(I,J)
                 A(I,J)=TEMP
               } // 22
               ENDIF
            } else {
              IP=-IPIV(I)
              if (I .GT. 1) {
              DO 23 J= 1,I-1
                 TEMP=A(IP,J)
                 A(IP,J)=A(I+1,J)
                 A(I+1,J)=TEMP
              } // 23
              ENDIF
              I=I+1
           ENDIF
           I=I+1
        END DO
         } else {

       // Revert A (A is lower)


         // Revert PERMUTATIONS

            I=N
            DO WHILE ( I .GE. 1 )
               if ( IPIV(I) .GT. 0 ) {
                  IP=IPIV(I)
                  if (I .GT. 1) {
                     DO J= 1,I-1
                        TEMP=A(I,J)
                        A(I,J)=A(IP,J)
                        A(IP,J)=TEMP
                     END DO
                  ENDIF
               } else {
                  IP=-IPIV(I)
                  I=I-1
                  if (I .GT. 1) {
                     DO J= 1,I-1
                        TEMP=A(I+1,J)
                        A(I+1,J)=A(IP,J)
                        A(IP,J)=TEMP
                     END DO
                  ENDIF
               ENDIF
               I=I-1
            END DO

         // Revert VALUE

            I=1
            DO WHILE ( I .LE. N-1 )
               if ( IPIV(I) .LT. 0 ) {
                  A(I+1,I)=E(I)
                  I=I+1
               ENDIF
               I=I+1
            END DO
         }
      }

      RETURN

      // End of CSYCONV

      }
