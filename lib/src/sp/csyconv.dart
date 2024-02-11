      void csyconv(final int UPLO, final int WAY, final int N, final Matrix<double> A, final int LDA, final Array<int> IPIV, final int E, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO, WAY;
      int                INFO, LDA, N;
      int                IPIV( * );
      Complex            A( LDA, * ), E( * );
      // ..

      Complex            ZERO;
      const              ZERO = (0.0,0.0) ;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame

      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // .. Local Scalars ..
      bool               UPPER, CONVERT;
      int                I, IP, J;
      Complex            TEMP;

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      CONVERT = lsame( WAY, 'C' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !CONVERT && !lsame( WAY, 'R' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;

      }
      if ( INFO != 0 ) {
         xerbla('CSYCONV', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( UPPER ) {

       // A is UPPER

       // Convert A (A is upper)

         // Convert VALUE

         if ( CONVERT ) {
            I=N;
            E(1)=ZERO;
            while (I > 1) {
               if ( IPIV(I) < 0 ) {
                  E(I)=A(I-1,I);
                  E(I-1)=ZERO;
                  A(I-1,I)=ZERO;
                  I=I-1;
               } else {
                  E(I)=ZERO;
               }
               I=I-1;
            }

         // Convert PERMUTATIONS

         I=N;
         while (I >= 1) {
            if ( IPIV(I) > 0) {
               IP=IPIV(I);
               if ( I < N) {
                  for (J = I+1; J <= N; J++) { // 12
                    TEMP=A(IP,J);
                    A(IP,J)=A(I,J);
                    A(I,J)=TEMP;
               } // 12
               }
            } else {
              IP=-IPIV(I);
               if ( I < N) {
             for (J = I+1; J <= N; J++) { // 13
                 TEMP=A(IP,J);
                 A(IP,J)=A(I-1,J);
                 A(I-1,J)=TEMP;
               } // 13
                }
                I=I-1;
           }
           I=I-1;
        }

         } else {

       // Revert A (A is upper)


         // Revert PERMUTATIONS

            I=1;
            while (I <= N) {
               if ( IPIV(I) > 0 ) {
                  IP=IPIV(I);
                  if ( I < N) {
                  for (J = I+1; J <= N; J++) {
                    TEMP=A(IP,J);
                    A(IP,J)=A(I,J);
                    A(I,J)=TEMP;
                  }
                  }
               } else {
                 IP=-IPIV(I);
                 I=I+1;
                 if ( I < N) {
                    for (J = I+1; J <= N; J++) {
                       TEMP=A(IP,J);
                       A(IP,J)=A(I-1,J);
                       A(I-1,J)=TEMP;
                    }
                 }
               }
               I=I+1;
            }

         // Revert VALUE

            I=N;
            while (I > 1) {
               if ( IPIV(I) < 0 ) {
                  A(I-1,I)=E(I);
                  I=I-1;
               }
               I=I-1;
            }
         }
      } else {

       // A is LOWER

         if ( CONVERT ) {

       // Convert A (A is lower)


         // Convert VALUE

            I=1;
            E(N)=ZERO;
            while (I <= N) {
               if ( I < N && IPIV(I) < 0 ) {
                  E(I)=A(I+1,I);
                  E(I+1)=ZERO;
                  A(I+1,I)=ZERO;
                  I=I+1;
               } else {
                  E(I)=ZERO;
               }
               I=I+1;
            }

         // Convert PERMUTATIONS

         I=1;
         while (I <= N) {
            if ( IPIV(I) > 0 ) {
               IP=IPIV(I);
               if (I > 1) {
               for (J = 1; J <= I-1; J++) { // 22
                 TEMP=A(IP,J);
                 A(IP,J)=A(I,J);
                 A(I,J)=TEMP;
               } // 22
               }
            } else {
              IP=-IPIV(I);
              if (I > 1) {
              for (J = 1; J <= I-1; J++) { // 23
                 TEMP=A(IP,J);
                 A(IP,J)=A(I+1,J);
                 A(I+1,J)=TEMP;
              } // 23
              }
              I=I+1;
           }
           I=I+1;
        }
         } else {

       // Revert A (A is lower)


         // Revert PERMUTATIONS

            I=N;
            while (I >= 1) {
               if ( IPIV(I) > 0 ) {
                  IP=IPIV(I);
                  if (I > 1) {
                     for (J = 1; J <= I-1; J++) {
                        TEMP=A(I,J);
                        A(I,J)=A(IP,J);
                        A(IP,J)=TEMP;
                     }
                  }
               } else {
                  IP=-IPIV(I);
                  I=I-1;
                  if (I > 1) {
                     for (J = 1; J <= I-1; J++) {
                        TEMP=A(I+1,J);
                        A(I+1,J)=A(IP,J);
                        A(IP,J)=TEMP;
                     }
                  }
               }
               I=I-1;
            }

         // Revert VALUE

            I=1;
            while (I <= N-1) {
               if ( IPIV(I) < 0 ) {
                  A(I+1,I)=E(I);
                  I=I+1;
               }
               I=I+1;
            }
         }
      }

      }
