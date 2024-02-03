      void slags2(UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, SNV, CSQ, SNQ ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               UPPER;
      REAL               A1, A2, A3, B1, B2, B3, CSQ, CSU, CSV, SNQ, SNU, SNV;
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      REAL               A, AUA11, AUA12, AUA21, AUA22, AVB11, AVB12, AVB21, AVB22, CSL, CSR, D, S1, S2, SNL, SNR, UA11R, UA22R, VB11R, VB22R, B, C, R, UA11, UA12, UA21, UA22, VB11, VB12, VB21, VB22;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARTG, SLASV2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      if ( UPPER ) {

         // Input matrices A and B are upper triangular matrices

         // Form matrix C = A*adj(B) = ( a b )
                                    // ( 0 d )

         A = A1*B3;
         D = A3*B1;
         B = A2*B1 - A1*B2;

         // The SVD of real 2-by-2 triangular C

          // ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 )
          // ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T )

         slasv2(A, B, D, S1, S2, SNR, CSR, SNL, CSL );

         if ( ABS( CSL ) >= ABS( SNL ) || ABS( CSR ) >= ABS( SNR ) ) {

            // Compute the (1,1) and (1,2) elements of U**T *A and V**T *B,
            // and (1,2) element of |U|**T *|A| and |V|**T *|B|.

            UA11R = CSL*A1;
            UA12 = CSL*A2 + SNL*A3;

            VB11R = CSR*B1;
            VB12 = CSR*B2 + SNR*B3;

            AUA12 = ABS( CSL )*ABS( A2 ) + ABS( SNL )*ABS( A3 );
            AVB12 = ABS( CSR )*ABS( B2 ) + ABS( SNR )*ABS( B3 );

            // zero (1,2) elements of U**T *A and V**T *B

            if ( ( ABS( UA11R )+ABS( UA12 ) ) != ZERO ) {
               if ( AUA12 / ( ABS( UA11R )+ABS( UA12 ) ) <= AVB12 / ( ABS( VB11R )+ABS( VB12 ) ) ) {
                  slartg(-UA11R, UA12, CSQ, SNQ, R );
               } else {
                  slartg(-VB11R, VB12, CSQ, SNQ, R );
               }
            } else {
               slartg(-VB11R, VB12, CSQ, SNQ, R );
            }

            CSU = CSL;
            SNU = -SNL;
            CSV = CSR;
            SNV = -SNR;

         } else {

            // Compute the (2,1) and (2,2) elements of U**T *A and V**T *B,
            // and (2,2) element of |U|**T *|A| and |V|**T *|B|.

            UA21 = -SNL*A1;
            UA22 = -SNL*A2 + CSL*A3;

            VB21 = -SNR*B1;
            VB22 = -SNR*B2 + CSR*B3;

            AUA22 = ABS( SNL )*ABS( A2 ) + ABS( CSL )*ABS( A3 );
            AVB22 = ABS( SNR )*ABS( B2 ) + ABS( CSR )*ABS( B3 );

            // zero (2,2) elements of U**T*A and V**T*B, and then swap.

            if ( ( ABS( UA21 )+ABS( UA22 ) ) != ZERO ) {
               if ( AUA22 / ( ABS( UA21 )+ABS( UA22 ) ) <= AVB22 / ( ABS( VB21 )+ABS( VB22 ) ) ) {
                  slartg(-UA21, UA22, CSQ, SNQ, R );
               } else {
                  slartg(-VB21, VB22, CSQ, SNQ, R );
               }
            } else {
               slartg(-VB21, VB22, CSQ, SNQ, R );
            }

            CSU = SNL;
            SNU = CSL;
            CSV = SNR;
            SNV = CSR;

         }

      } else {

         // Input matrices A and B are lower triangular matrices

         // Form matrix C = A*adj(B) = ( a 0 )
                                    // ( c d )

         A = A1*B3;
         D = A3*B1;
         C = A2*B3 - A3*B2;

         // The SVD of real 2-by-2 triangular C

          // ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 )
          // ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T )

         slasv2(A, C, D, S1, S2, SNR, CSR, SNL, CSL );

         if ( ABS( CSR ) >= ABS( SNR ) || ABS( CSL ) >= ABS( SNL ) ) {

            // Compute the (2,1) and (2,2) elements of U**T *A and V**T *B,
            // and (2,1) element of |U|**T *|A| and |V|**T *|B|.

            UA21 = -SNR*A1 + CSR*A2;
            UA22R = CSR*A3;

            VB21 = -SNL*B1 + CSL*B2;
            VB22R = CSL*B3;

            AUA21 = ABS( SNR )*ABS( A1 ) + ABS( CSR )*ABS( A2 );
            AVB21 = ABS( SNL )*ABS( B1 ) + ABS( CSL )*ABS( B2 );

            // zero (2,1) elements of U**T *A and V**T *B.

            if ( ( ABS( UA21 )+ABS( UA22R ) ) != ZERO ) {
               if ( AUA21 / ( ABS( UA21 )+ABS( UA22R ) ) <= AVB21 / ( ABS( VB21 )+ABS( VB22R ) ) ) {
                  slartg(UA22R, UA21, CSQ, SNQ, R );
               } else {
                  slartg(VB22R, VB21, CSQ, SNQ, R );
               }
            } else {
               slartg(VB22R, VB21, CSQ, SNQ, R );
            }

            CSU = CSR;
            SNU = -SNR;
            CSV = CSL;
            SNV = -SNL;

         } else {

            // Compute the (1,1) and (1,2) elements of U**T *A and V**T *B,
            // and (1,1) element of |U|**T *|A| and |V|**T *|B|.

            UA11 = CSR*A1 + SNR*A2;
            UA12 = SNR*A3;

            VB11 = CSL*B1 + SNL*B2;
            VB12 = SNL*B3;

            AUA11 = ABS( CSR )*ABS( A1 ) + ABS( SNR )*ABS( A2 );
            AVB11 = ABS( CSL )*ABS( B1 ) + ABS( SNL )*ABS( B2 );

            // zero (1,1) elements of U**T*A and V**T*B, and then swap.

            if ( ( ABS( UA11 )+ABS( UA12 ) ) != ZERO ) {
               if ( AUA11 / ( ABS( UA11 )+ABS( UA12 ) ) <= AVB11 / ( ABS( VB11 )+ABS( VB12 ) ) ) {
                  slartg(UA12, UA11, CSQ, SNQ, R );
               } else {
                  slartg(VB12, VB11, CSQ, SNQ, R );
               }
            } else {
               slartg(VB12, VB11, CSQ, SNQ, R );
            }

            CSU = SNR;
            SNU = CSR;
            CSV = SNL;
            SNV = CSL;

         }

      }

      return;
      }
