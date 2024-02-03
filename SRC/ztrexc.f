      SUBROUTINE ZTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPQ;
      int                IFST, ILST, INFO, LDQ, LDT, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         Q( LDQ, * ), T( LDT, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               WANTQ;
      int                K, M1, M2, M3;
      double             CS;
      COMPLEX*16         SN, T11, T22, TEMP
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLARTG, ZROT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG, MAX
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters.

      INFO = 0
      WANTQ = LSAME( COMPQ, 'V' )
      IF( .NOT.LSAME( COMPQ, 'N' ) .AND. .NOT.WANTQ ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.MAX( 1, N ) ) ) THEN
         INFO = -6
      ELSE IF(( IFST.LT.1 .OR. IFST.GT.N ).AND.( N.GT.0 )) THEN
         INFO = -7
      ELSE IF(( ILST.LT.1 .OR. ILST.GT.N ).AND.( N.GT.0 )) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTREXC', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.LE.1 .OR. IFST.EQ.ILST ) RETURN

      IF( IFST.LT.ILST ) THEN

         // Move the IFST-th diagonal element forward down the diagonal.

         M1 = 0
         M2 = -1
         M3 = 1
      ELSE

         // Move the IFST-th diagonal element backward up the diagonal.

         M1 = -1
         M2 = 0
         M3 = -1
      END IF

      DO 10 K = IFST + M1, ILST + M2, M3

         // Interchange the k-th and (k+1)-th diagonal elements.

         T11 = T( K, K )
         T22 = T( K+1, K+1 )

         // Determine the transformation to perform the interchange.

         CALL ZLARTG( T( K, K+1 ), T22-T11, CS, SN, TEMP )

         // Apply transformation to the matrix T.

         IF( K+2.LE.N ) CALL ZROT( N-K-1, T( K, K+2 ), LDT, T( K+1, K+2 ), LDT, CS, SN )
         CALL ZROT( K-1, T( 1, K ), 1, T( 1, K+1 ), 1, CS, DCONJG( SN ) )

         T( K, K ) = T22
         T( K+1, K+1 ) = T11

         IF( WANTQ ) THEN

            // Accumulate transformation in the matrix Q.

            CALL ZROT( N, Q( 1, K ), 1, Q( 1, K+1 ), 1, CS, DCONJG( SN ) )
         END IF

   10 CONTINUE

      RETURN

      // End of ZTREXC

      }
