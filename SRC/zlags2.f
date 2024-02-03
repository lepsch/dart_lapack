      SUBROUTINE ZLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, SNV, CSQ, SNQ )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               UPPER;
      double             A1, A3, B1, B3, CSQ, CSU, CSV;
      COMPLEX*16         A2, B2, SNQ, SNU, SNV
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      double             A, AUA11, AUA12, AUA21, AUA22, AVB12, AVB11, AVB21, AVB22, CSL, CSR, D, FB, FC, S1, S2, SNL, SNR, UA11R, UA22R, VB11R, VB22R;
      COMPLEX*16         B, C, D1, R, T, UA11, UA12, UA21, UA22, VB11, VB12, VB21, VB22
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASV2, ZLARTG
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, DIMAG
      // ..
      // .. Statement Functions ..
      double             ABS1;
      // ..
      // .. Statement Function definitions ..
      ABS1( T ) = ABS( DBLE( T ) ) + ABS( DIMAG( T ) )
      // ..
      // .. Executable Statements ..

      if ( UPPER ) {

         // Input matrices A and B are upper triangular matrices

         // Form matrix C = A*adj(B) = ( a b )
                                    // ( 0 d )

         A = A1*B3
         D = A3*B1
         B = A2*B1 - A1*B2
         FB = ABS( B )

         // Transform complex 2-by-2 matrix C to real matrix by unitary
         // diagonal matrix diag(1,D1).

         D1 = ONE
         IF( FB.NE.ZERO ) D1 = B / FB

         // The SVD of real 2 by 2 triangular C

          // ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 )
          // ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T )

         CALL DLASV2( A, FB, D, S1, S2, SNR, CSR, SNL, CSL )

         if ( ABS( CSL ).GE.ABS( SNL ) .OR. ABS( CSR ).GE.ABS( SNR ) ) {

            // Compute the (1,1) and (1,2) elements of U**H *A and V**H *B,
            // and (1,2) element of |U|**H *|A| and |V|**H *|B|.

            UA11R = CSL*A1
            UA12 = CSL*A2 + D1*SNL*A3

            VB11R = CSR*B1
            VB12 = CSR*B2 + D1*SNR*B3

            AUA12 = ABS( CSL )*ABS1( A2 ) + ABS( SNL )*ABS( A3 )
            AVB12 = ABS( CSR )*ABS1( B2 ) + ABS( SNR )*ABS( B3 )

            // zero (1,2) elements of U**H *A and V**H *B

            if ( ( ABS( UA11R )+ABS1( UA12 ) ).EQ.ZERO ) {
               CALL ZLARTG( -DCMPLX( VB11R ), DCONJG( VB12 ), CSQ, SNQ, R )
            } else if ( ( ABS( VB11R )+ABS1( VB12 ) ).EQ.ZERO ) {
               CALL ZLARTG( -DCMPLX( UA11R ), DCONJG( UA12 ), CSQ, SNQ, R )             ELSE IF( AUA12 / ( ABS( UA11R )+ABS1( UA12 ) ).LE.AVB12 / ( ABS( VB11R )+ABS1( VB12 ) ) ) THEN                CALL ZLARTG( -DCMPLX( UA11R ), DCONJG( UA12 ), CSQ, SNQ, R )
            } else {
               CALL ZLARTG( -DCMPLX( VB11R ), DCONJG( VB12 ), CSQ, SNQ, R )
            }

            CSU = CSL
            SNU = -D1*SNL
            CSV = CSR
            SNV = -D1*SNR

         } else {

            // Compute the (2,1) and (2,2) elements of U**H *A and V**H *B,
            // and (2,2) element of |U|**H *|A| and |V|**H *|B|.

            UA21 = -DCONJG( D1 )*SNL*A1
            UA22 = -DCONJG( D1 )*SNL*A2 + CSL*A3

            VB21 = -DCONJG( D1 )*SNR*B1
            VB22 = -DCONJG( D1 )*SNR*B2 + CSR*B3

            AUA22 = ABS( SNL )*ABS1( A2 ) + ABS( CSL )*ABS( A3 )
            AVB22 = ABS( SNR )*ABS1( B2 ) + ABS( CSR )*ABS( B3 )

            // zero (2,2) elements of U**H *A and V**H *B, and then swap.

            if ( ( ABS1( UA21 )+ABS1( UA22 ) ).EQ.ZERO ) {
               CALL ZLARTG( -DCONJG( VB21 ), DCONJG( VB22 ), CSQ, SNQ, R )
            } else if ( ( ABS1( VB21 )+ABS( VB22 ) ).EQ.ZERO ) {
               CALL ZLARTG( -DCONJG( UA21 ), DCONJG( UA22 ), CSQ, SNQ, R )             ELSE IF( AUA22 / ( ABS1( UA21 )+ABS1( UA22 ) ).LE.AVB22 / ( ABS1( VB21 )+ABS1( VB22 ) ) ) THEN                CALL ZLARTG( -DCONJG( UA21 ), DCONJG( UA22 ), CSQ, SNQ, R )
            } else {
               CALL ZLARTG( -DCONJG( VB21 ), DCONJG( VB22 ), CSQ, SNQ, R )
            }

            CSU = SNL
            SNU = D1*CSL
            CSV = SNR
            SNV = D1*CSR

         }

      } else {

         // Input matrices A and B are lower triangular matrices

         // Form matrix C = A*adj(B) = ( a 0 )
                                    // ( c d )

         A = A1*B3
         D = A3*B1
         C = A2*B3 - A3*B2
         FC = ABS( C )

         // Transform complex 2-by-2 matrix C to real matrix by unitary
         // diagonal matrix diag(d1,1).

         D1 = ONE
         IF( FC.NE.ZERO ) D1 = C / FC

         // The SVD of real 2 by 2 triangular C

          // ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 )
          // ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T )

         CALL DLASV2( A, FC, D, S1, S2, SNR, CSR, SNL, CSL )

         if ( ABS( CSR ).GE.ABS( SNR ) .OR. ABS( CSL ).GE.ABS( SNL ) ) {

            // Compute the (2,1) and (2,2) elements of U**H *A and V**H *B,
            // and (2,1) element of |U|**H *|A| and |V|**H *|B|.

            UA21 = -D1*SNR*A1 + CSR*A2
            UA22R = CSR*A3

            VB21 = -D1*SNL*B1 + CSL*B2
            VB22R = CSL*B3

            AUA21 = ABS( SNR )*ABS( A1 ) + ABS( CSR )*ABS1( A2 )
            AVB21 = ABS( SNL )*ABS( B1 ) + ABS( CSL )*ABS1( B2 )

            // zero (2,1) elements of U**H *A and V**H *B.

            if ( ( ABS1( UA21 )+ABS( UA22R ) ).EQ.ZERO ) {
               CALL ZLARTG( DCMPLX( VB22R ), VB21, CSQ, SNQ, R )
            } else if ( ( ABS1( VB21 )+ABS( VB22R ) ).EQ.ZERO ) {
               CALL ZLARTG( DCMPLX( UA22R ), UA21, CSQ, SNQ, R )
            } else if ( AUA21 / ( ABS1( UA21 )+ABS( UA22R ) ).LE.AVB21 / ( ABS1( VB21 )+ABS( VB22R ) ) ) {
               CALL ZLARTG( DCMPLX( UA22R ), UA21, CSQ, SNQ, R )
            } else {
               CALL ZLARTG( DCMPLX( VB22R ), VB21, CSQ, SNQ, R )
            }

            CSU = CSR
            SNU = -DCONJG( D1 )*SNR
            CSV = CSL
            SNV = -DCONJG( D1 )*SNL

         } else {

            // Compute the (1,1) and (1,2) elements of U**H *A and V**H *B,
            // and (1,1) element of |U|**H *|A| and |V|**H *|B|.

            UA11 = CSR*A1 + DCONJG( D1 )*SNR*A2
            UA12 = DCONJG( D1 )*SNR*A3

            VB11 = CSL*B1 + DCONJG( D1 )*SNL*B2
            VB12 = DCONJG( D1 )*SNL*B3

            AUA11 = ABS( CSR )*ABS( A1 ) + ABS( SNR )*ABS1( A2 )
            AVB11 = ABS( CSL )*ABS( B1 ) + ABS( SNL )*ABS1( B2 )

            // zero (1,1) elements of U**H *A and V**H *B, and then swap.

            if ( ( ABS1( UA11 )+ABS1( UA12 ) ).EQ.ZERO ) {
               CALL ZLARTG( VB12, VB11, CSQ, SNQ, R )
            } else if ( ( ABS1( VB11 )+ABS1( VB12 ) ).EQ.ZERO ) {
               CALL ZLARTG( UA12, UA11, CSQ, SNQ, R )
            } else if ( AUA11 / ( ABS1( UA11 )+ABS1( UA12 ) ).LE.AVB11 / ( ABS1( VB11 )+ABS1( VB12 ) ) ) {
               CALL ZLARTG( UA12, UA11, CSQ, SNQ, R )
            } else {
               CALL ZLARTG( VB12, VB11, CSQ, SNQ, R )
            }

            CSU = SNR
            SNU = DCONJG( D1 )*CSR
            CSV = SNL
            SNV = DCONJG( D1 )*CSL

         }

      }

      RETURN

      // End of ZLAGS2

      }
