      void clags2(UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, SNV, CSQ, SNQ ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               UPPER;
      double               A1, A3, B1, B3, CSQ, CSU, CSV;
      Complex            A2, B2, SNQ, SNU, SNV;
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      double               A, AUA11, AUA12, AUA21, AUA22, AVB11, AVB12, AVB21, AVB22, CSL, CSR, D, FB, FC, S1, S2, SNL, SNR, UA11R, UA22R, VB11R, VB22R;
      Complex            B, C, D1, R, T, UA11, UA12, UA21, UA22, VB11, VB12, VB21, VB22;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARTG, SLASV2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, CONJG, REAL
      // ..
      // .. Statement Functions ..
      double               ABS1;
      // ..
      // .. Statement Function definitions ..
      ABS1[T] = ( double( T ) ).abs() + ( AIMAG( T ) ).abs();
      // ..
      // .. Executable Statements ..

      if ( UPPER ) {

         // Input matrices A and B are upper triangular matrices

         // Form matrix C = A*adj(B) = ( a b )
                                    // ( 0 d )

         A = A1*B3;
         D = A3*B1;
         B = A2*B1 - A1*B2;
         FB = ( B ).abs();

         // Transform complex 2-by-2 matrix C to real matrix by unitary
         // diagonal matrix diag(1,D1).

         D1 = ONE;
         if (FB != ZERO) D1 = B / FB;

         // The SVD of real 2 by 2 triangular C

          // ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 )
          // ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T )

         slasv2(A, FB, D, S1, S2, SNR, CSR, SNL, CSL );

         if ( ( CSL ).abs() >= ( SNL ).abs() || ( CSR ).abs() >= ( SNR ).abs() ) {

            // Compute the (1,1) and (1,2) elements of U**H *A and V**H *B,
            // and (1,2) element of |U|**H *|A| and |V|**H *|B|.

            UA11R = CSL*A1;
            UA12 = CSL*A2 + D1*SNL*A3;

            VB11R = CSR*B1;
            VB12 = CSR*B2 + D1*SNR*B3;

            AUA12 = ( CSL ).abs()*ABS1( A2 ) + ( SNL ).abs()*( A3 ).abs();
            AVB12 = ( CSR ).abs()*ABS1( B2 ) + ( SNR ).abs()*( B3 ).abs();

            // zero (1,2) elements of U**H *A and V**H *B

            if ( ( ( UA11R ).abs()+ABS1( UA12 ) ) == ZERO ) {
               clartg(-CMPLX( VB11R ), CONJG( VB12 ), CSQ, SNQ, R );
            } else if ( ( ( VB11R ).abs()+ABS1( VB12 ) ) == ZERO ) {
               clartg(-CMPLX( UA11R ), CONJG( UA12 ), CSQ, SNQ, R )             ELSE IF( AUA12 / ( ( UA11R ).abs()+ABS1( UA12 ) ) <= AVB12 / ( ( VB11R ).abs()+ABS1( VB12 ) ) ) THEN;
               clartg(-CMPLX( UA11R ), CONJG( UA12 ), CSQ, SNQ, R );
            } else {
               clartg(-CMPLX( VB11R ), CONJG( VB12 ), CSQ, SNQ, R );
            }

            CSU = CSL;
            SNU = -D1*SNL;
            CSV = CSR;
            SNV = -D1*SNR;

         } else {

            // Compute the (2,1) and (2,2) elements of U**H *A and V**H *B,
            // and (2,2) element of |U|**H *|A| and |V|**H *|B|.

            UA21 = -CONJG( D1 )*SNL*A1;
            UA22 = -CONJG( D1 )*SNL*A2 + CSL*A3;

            VB21 = -CONJG( D1 )*SNR*B1;
            VB22 = -CONJG( D1 )*SNR*B2 + CSR*B3;

            AUA22 = ( SNL ).abs()*ABS1( A2 ) + ( CSL ).abs()*( A3 ).abs();
            AVB22 = ( SNR ).abs()*ABS1( B2 ) + ( CSR ).abs()*( B3 ).abs();

            // zero (2,2) elements of U**H *A and V**H *B, and then swap.

            if ( ( ABS1( UA21 )+ABS1( UA22 ) ) == ZERO ) {
               clartg(-CONJG( VB21 ), CONJG( VB22 ), CSQ, SNQ, R );
            } else if ( ( ABS1( VB21 )+( VB22 ).abs() ) == ZERO ) {
               clartg(-CONJG( UA21 ), CONJG( UA22 ), CSQ, SNQ, R );
            } else if ( AUA22 / ( ABS1( UA21 )+ABS1( UA22 ) ) <= AVB22 / ( ABS1( VB21 )+ABS1( VB22 ) ) ) {
               clartg(-CONJG( UA21 ), CONJG( UA22 ), CSQ, SNQ, R );
            } else {
               clartg(-CONJG( VB21 ), CONJG( VB22 ), CSQ, SNQ, R );
            }

            CSU = SNL;
            SNU = D1*CSL;
            CSV = SNR;
            SNV = D1*CSR;

         }

      } else {

         // Input matrices A and B are lower triangular matrices

         // Form matrix C = A*adj(B) = ( a 0 )
                                    // ( c d )

         A = A1*B3;
         D = A3*B1;
         C = A2*B3 - A3*B2;
         FC = ( C ).abs();

         // Transform complex 2-by-2 matrix C to real matrix by unitary
         // diagonal matrix diag(d1,1).

         D1 = ONE;
         if (FC != ZERO) D1 = C / FC;

         // The SVD of real 2 by 2 triangular C

          // ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 )
          // ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T )

         slasv2(A, FC, D, S1, S2, SNR, CSR, SNL, CSL );

         if ( ( CSR ).abs() >= ( SNR ).abs() || ( CSL ).abs() >= ( SNL ).abs() ) {

            // Compute the (2,1) and (2,2) elements of U**H *A and V**H *B,
            // and (2,1) element of |U|**H *|A| and |V|**H *|B|.

            UA21 = -D1*SNR*A1 + CSR*A2;
            UA22R = CSR*A3;

            VB21 = -D1*SNL*B1 + CSL*B2;
            VB22R = CSL*B3;

            AUA21 = ( SNR ).abs()*( A1 ).abs() + ( CSR ).abs()*ABS1( A2 );
            AVB21 = ( SNL ).abs()*( B1 ).abs() + ( CSL ).abs()*ABS1( B2 );

            // zero (2,1) elements of U**H *A and V**H *B.

            if ( ( ABS1( UA21 )+( UA22R ).abs() ) == ZERO ) {
               clartg(CMPLX( VB22R ), VB21, CSQ, SNQ, R );
            } else if ( ( ABS1( VB21 )+( VB22R ).abs() ) == ZERO ) {
               clartg(CMPLX( UA22R ), UA21, CSQ, SNQ, R );
            } else if ( AUA21 / ( ABS1( UA21 )+( UA22R ).abs() ) <= AVB21 / ( ABS1( VB21 )+( VB22R ).abs() ) ) {
               clartg(CMPLX( UA22R ), UA21, CSQ, SNQ, R );
            } else {
               clartg(CMPLX( VB22R ), VB21, CSQ, SNQ, R );
            }

            CSU = CSR;
            SNU = -CONJG( D1 )*SNR;
            CSV = CSL;
            SNV = -CONJG( D1 )*SNL;

         } else {

            // Compute the (1,1) and (1,2) elements of U**H *A and V**H *B,
            // and (1,1) element of |U|**H *|A| and |V|**H *|B|.

            UA11 = CSR*A1 + CONJG( D1 )*SNR*A2;
            UA12 = CONJG( D1 )*SNR*A3;

            VB11 = CSL*B1 + CONJG( D1 )*SNL*B2;
            VB12 = CONJG( D1 )*SNL*B3;

            AUA11 = ( CSR ).abs()*( A1 ).abs() + ( SNR ).abs()*ABS1( A2 );
            AVB11 = ( CSL ).abs()*( B1 ).abs() + ( SNL ).abs()*ABS1( B2 );

            // zero (1,1) elements of U**H *A and V**H *B, and then swap.

            if ( ( ABS1( UA11 )+ABS1( UA12 ) ) == ZERO ) {
               clartg(VB12, VB11, CSQ, SNQ, R );
            } else if ( ( ABS1( VB11 )+ABS1( VB12 ) ) == ZERO ) {
               clartg(UA12, UA11, CSQ, SNQ, R );
            } else if ( AUA11 / ( ABS1( UA11 )+ABS1( UA12 ) ) <= AVB11 / ( ABS1( VB11 )+ABS1( VB12 ) ) ) {
               clartg(UA12, UA11, CSQ, SNQ, R );
            } else {
               clartg(VB12, VB11, CSQ, SNQ, R );
            }

            CSU = SNR;
            SNU = CONJG( D1 )*CSR;
            CSV = SNL;
            SNV = CONJG( D1 )*CSL;

         }

      }

      return;
      }