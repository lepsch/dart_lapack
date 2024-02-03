      SUBROUTINE CLAHILB( N, NRHS, A, LDA, X, LDX, B, LDB, WORK, INFO, PATH)

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     N, NRHS, LDA, LDX, LDB, INFO;
      // .. Array Arguments ..
      REAL WORK(N)
      COMPLEX A(LDA,N), X(LDX, NRHS), B(LDB, NRHS)
      String             PATH;
      // ..

*  =====================================================================
      // .. Local Scalars ..
      int     TM, TI, R;
      int     M;
      int     I, J;
      COMPLEX TMP
      String      C2;
      // ..
      // .. Parameters ..
      // NMAX_EXACT   the largest dimension where the generated data is
                   // exact.
      // NMAX_APPROX  the largest dimension where the generated data has
                   // a small componentwise relative error.
      // ??? complex uses how many bits ???
      int     NMAX_EXACT, NMAX_APPROX, SIZE_D;
      const     NMAX_EXACT = 6, NMAX_APPROX = 11, SIZE_D = 8;

      // d's are generated from random permutation of those eight elements.
      COMPLEX D1(8), D2(8), INVD1(8), INVD2(8)
      DATA D1 /(-1,0),(0,1),(-1,-1),(0,-1),(1,0),(-1,1),(1,1),(1,-1)/
      DATA D2 /(-1,0),(0,-1),(-1,1),(0,1),(1,0),(-1,-1),(1,-1),(1,1)/
       DATA INVD1 /(-1,0),(0,-1),(-.5,.5),(0,1),(1,0), (-.5,-.5),(.5,-.5),(.5,.5)/       DATA INVD2 /(-1,0),(0,1),(-.5,-.5),(0,-1),(1,0), (-.5,.5),(.5,.5),(.5,-.5)/
      // ..
      // .. External Functions
      // EXTERNAL CLASET, LSAMEN
      // INTRINSIC REAL
      bool    LSAMEN;
      // ..
      // .. Executable Statements ..
      C2 = PATH( 2: 3 )

      // Test the input arguments

      INFO = 0
      if (N .LT. 0 .OR. N .GT. NMAX_APPROX) {
         INFO = -1
      } else if (NRHS .LT. 0) {
         INFO = -2
      } else if (LDA .LT. N) {
         INFO = -4
      } else if (LDX .LT. N) {
         INFO = -6
      } else if (LDB .LT. N) {
         INFO = -8
      }
      if (INFO .LT. 0) {
         xerbla('CLAHILB', -INFO);
         RETURN
      }
      if (N .GT. NMAX_EXACT) {
         INFO = 1
      }

      // Compute M = the LCM of the integers [1, 2*N-1].  The largest
      // reasonable N is small enough that integers suffice (up to N = 11).
      M = 1
      for (I = 2; I <= (2*N-1); I++) {
         TM = M
         TI = I
         R = MOD(TM, TI)
         DO WHILE (R .NE. 0)
            TM = TI
            TI = R
            R = MOD(TM, TI)
         END DO
         M = (M / TI) * I
      END DO

      // Generate the scaled Hilbert matrix in A
      // If we are testing SY routines, take
           // D1_i = D2_i, else, D1_i = D2_i*
      if ( LSAMEN( 2, C2, 'SY' ) ) {
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= N; I++) {
               A(I, J) = D1(MOD(J,SIZE_D)+1) * (REAL(M) / (I + J - 1)) * D1(MOD(I,SIZE_D)+1)
            END DO
         END DO
      } else {
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= N; I++) {
               A(I, J) = D1(MOD(J,SIZE_D)+1) * (REAL(M) / (I + J - 1)) * D2(MOD(I,SIZE_D)+1)
            END DO
         END DO
      }

      // Generate matrix B as simply the first NRHS columns of M * the
      // identity.
      TMP = REAL(M)
      claset('Full', N, NRHS, (0.0,0.0), TMP, B, LDB);

      // Generate the true solutions in X.  Because B = the first NRHS
      // columns of M*I, the true solutions are just the first NRHS columns
      // of the inverse Hilbert matrix.
      WORK(1) = N
      for (J = 2; J <= N; J++) {
         WORK(J) = (  ( (WORK(J-1)/(J-1)) * (J-1 - N) ) /(J-1)  ) * (N +J -1)
      END DO

      // If we are testing SY routines,
             // take D1_i = D2_i, else, D1_i = D2_i*
      if ( LSAMEN( 2, C2, 'SY' ) ) {
         for (J = 1; J <= NRHS; J++) {
            for (I = 1; I <= N; I++) {
               X(I, J) = INVD1(MOD(J,SIZE_D)+1) * ((WORK(I)*WORK(J)) / (I + J - 1)) * INVD1(MOD(I,SIZE_D)+1)
            END DO
         END DO
      } else {
         for (J = 1; J <= NRHS; J++) {
            for (I = 1; I <= N; I++) {
               X(I, J) = INVD2(MOD(J,SIZE_D)+1) * ((WORK(I)*WORK(J)) / (I + J - 1)) * INVD1(MOD(I,SIZE_D)+1)
            END DO
         END DO
      }
      }
