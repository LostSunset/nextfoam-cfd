/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | baramPODreconstruct - folder-scan + interpolation
   \\    /   O peration     | listField.dat: <field> <type> <regionName|->
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "IOmanip.H"
#include "IFstream.H"
#include "OFstream.H"
#include "processorFvPatch.H"
#include "fvPatchField.H"
#include "coupledFvPatchField.H"

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif
#include <Eigen/Dense>
#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

using namespace Foam;
using namespace Eigen;

namespace
{

struct FieldSpec
{
    word name;
    word type;   // "scalar" or "vector"
    word region; // regionName or "-"
};

static inline bool isVector(const word& t)
{
    return (t == "vector" || t == "Vector" || t == "VECTOR" || t == "v" || t == "V");
}

// listField.dat: <field> <type> <region|->
DynamicList<FieldSpec> readFieldList(const fileName& path)
{
    DynamicList<FieldSpec> out;
    IFstream is(path);
    if (!is.good())
    {
        FatalErrorInFunction << "Cannot open " << path << exit(FatalError);
    }

    string line;
    DynamicList<FieldSpec> acc;
    while (is.getLine(line))
    {
        auto trim = [](string& s){
            while (!s.empty() && (s.back()=='\r' || s.back()=='\n' || s.back()==' ' || s.back()=='\t')) s.pop_back();
            size_t p=0; while (p<s.size() && (s[p]==' '||s[p]=='\t')) ++p; if (p) s.erase(0,p);
        };
        trim(line);
        if (line.empty() || line[0]=='#') continue;

        std::istringstream ss(line);
        FieldSpec fs;
        if (!(ss >> fs.name >> fs.type >> fs.region))
        {
            FatalErrorInFunction
                << "Invalid line in " << path << ": \"" << line << "\"\n"
                << "Expected format: <fieldName> <type> <regionName|->" << exit(FatalError);
        }
        acc.append(fs);
    }
    out.transfer(acc);
    return out;
}

// ---------- Flexible ASCII matrix readers ----------
// Format A: first line "rows cols", then rows lines of numbers
// Format B: plain whitespace matrix (rows unknown, cols = first line cols)
MatrixXd readMatrixFlexible(const fileName& fn)
{
    // Try format A
    {
        IFstream is(fn);
        if (!is.good())
            FatalErrorInFunction << "Cannot read: " << fn << exit(FatalError);

        string first;
        if (is.getLine(first))
        {
            std::istringstream hdr(first);
            label r, c;
            if (hdr >> r >> c)
            {
                MatrixXd M(r, c);
                string line;
                for (label i=0;i<r;++i)
                {
                    if (!is.getLine(line))
                        FatalErrorInFunction << "Unexpected EOF in " << fn << " at row " << i << exit(FatalError);
                    std::istringstream ss(line);
                    for (label j=0;j<c;++j)
                    {
                        double v; if (!(ss >> v))
                            FatalErrorInFunction << "Bad value in " << fn << " at ("<<i<<","<<j<<")" << exit(FatalError);
                        M(i,j) = v;
                    }
                }
                return M;
            }
        }
    }

    // Fallback: format B (plain)
    {
        IFstream is(fn);
        if (!is.good())
            FatalErrorInFunction << "Cannot read: " << fn << exit(FatalError);

        string line;
        std::vector<std::vector<double>> rows;
        while (is.getLine(line))
        {
            // skip empty/comment
            bool allws=true;
            for (char ch: line){ if (!std::isspace(static_cast<unsigned char>(ch))) {allws=false; break;} }
            if (allws || (!line.empty() && line[0]=='#')) continue;

            std::istringstream ss(line);
            std::vector<double> row;
            double v;
            while (ss >> v) row.push_back(v);
            if (!row.empty()) rows.push_back(std::move(row));
        }
        if (rows.empty()) FatalErrorInFunction << "Empty matrix file: " << fn << exit(FatalError);

        const label c = static_cast<label>(rows.front().size());
        for (const auto& r : rows)
            if (static_cast<label>(r.size()) != c)
                FatalErrorInFunction << "Inconsistent columns in " << fn << exit(FatalError);

        const label r = static_cast<label>(rows.size());
        MatrixXd M(r, c);
        for (label i=0;i<r;++i)
            for (label j=0;j<c;++j)
                M(i,j) = rows[i][j];
        return M;
    }
}

VectorXd readVectorFlexible(const fileName& fn)
{
    MatrixXd M = readMatrixFlexible(fn);
    if (M.cols()==1) return M.col(0);
    if (M.rows()==1) return M.row(0).transpose();
    FatalErrorInFunction << "Vector expected in " << fn << " but got "
                         << M.rows() << "x" << M.cols() << exit(FatalError);
    return VectorXd(); // unreachable
}

static inline dimensionSet dimsOf(const word& fld, bool isVector)
{
    if (isVector)           return dimensionSet(0, 1, -1, 0, 0, 0, 0);       // U차원
    if (fld == "rho")       return dimensionSet(1, -3, 0, 0, 0, 0, 0);       // 밀도
    if (fld == "p")         return dimensionSet(1, -1, -2, 0, 0, 0, 0);      // 압력
    if (fld == "T")         return dimensionSet(0, 0, 0, 1, 0, 0, 0);        // 온도
    if (fld == "Umag")      return dimensionSet(0, 1, -1, 0, 0, 0, 0);       // 속도 크기
    return dimensionSet(0, 0, 0, 0, 0, 0, 0);                                 // 디폴트 무차원
}

static void updateOnlyProcessorBC(volScalarField& fld)
{
    auto& b = fld.boundaryFieldRef();
    forAll(b, patchI)
    {
        const fvPatch& p = b[patchI].patch();
        if (isA<processorFvPatch>(p))
        {
            b[patchI].updateCoeffs();
            if (isA<coupledFvPatchField<scalar>>(b[patchI]))
            {
                auto& c = refCast<coupledFvPatchField<scalar>>(b[patchI]);
                const auto ct = Foam::UPstream::commsTypes::blocking;
                c.initEvaluate(ct); c.evaluate(ct);
            }
        }
    }
}

static void updateOnlyProcessorBC(volVectorField& fld)
{
    auto& b = fld.boundaryFieldRef();
    forAll(b, patchI)
    {
        const fvPatch& p = b[patchI].patch();
        if (isA<processorFvPatch>(p))
        {
            b[patchI].updateCoeffs();
            if (isA<coupledFvPatchField<vector>>(b[patchI]))
            {
                auto& c = refCast<coupledFvPatchField<vector>>(b[patchI]);
                const auto ct = Foam::UPstream::commsTypes::blocking;
                c.initEvaluate(ct); c.evaluate(ct);
            }
        }
    }
}

// --- 폴더 스캔: constant/*/polyMesh/points
static void scanRegions(const Time& runTime, wordList& regions, bool& multiRegion)
{
    regions.clear();
    const fileName constDir = runTime.constant();

    // 1) 디렉터리 전용
    fileNameList dirNames = readDir(constDir, fileName::DIRECTORY);

    // 2) 폴백: 전체→isDir 필터
    if (dirNames.empty())
    {
        fileNameList all = readDir(constDir);
        forAll(all, i)
        {
            const fileName sub = constDir/all[i];
            if (isDir(sub)) dirNames.append(all[i]);
        }
    }

    // 3) polyMesh/points 있는 폴더만 채택
    forAll(dirNames, i)
    {
        const fileName pointsFile = constDir/dirNames[i]/"polyMesh"/"points";
        if (isFile(pointsFile)) regions.append(dirNames[i]);
    }

    multiRegion = !regions.empty();
}

static label nEntriesWithBoundary(const volScalarField& s)
{
    label n = s.internalField().size();
    const auto& bf = s.boundaryField();
    forAll(bf, pi)
        if (!isA<processorFvPatch>(bf[pi].patch()))
            n += bf[pi].size();
    return n;
}

static label nEntriesWithBoundary(const volVectorField& v)
{
    label n = v.internalField().size();
    const auto& bf = v.boundaryField();
    forAll(bf, pi)
        if (!isA<processorFvPatch>(bf[pi].patch()))
            n += bf[pi].size();
    return 3*n;
}

void writeScalarToField(fvMesh& mesh, const word& fld, const Eigen::VectorXd& y)
{
    auto makeNew = [&]()
    {
        return volScalarField
        (
            IOobject(fld, mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
            mesh,
            dimensionedScalar(fld, dimsOf(fld, /*isVector*/false), 0.0)
        );
    };

    volScalarField* sp = nullptr;
    {
        IOobject hdr(fld, "0", mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE);
        sp = hdr.typeHeaderOk<volScalarField>(true) ? new volScalarField(hdr, mesh)
                                                    : new volScalarField(makeNew());
    }
    auto& s = *sp;

    if (static_cast<label>(y.size()) != nEntriesWithBoundary(s))
        FatalErrorInFunction << "Size mismatch writing " << fld << exit(FatalError);

    label k = 0;
    forAll(s, i) s[i] = y(k++);                        // internal
    auto& bf = s.boundaryFieldRef();
    forAll(bf, pi)                                     // boundary (non-processor)
        if (!isA<processorFvPatch>(bf[pi].patch()))
            forAll(bf[pi], f) bf[pi][f] = y(k++);

    updateOnlyProcessorBC(s);                          // processor 동기화만
    s.write();
    delete sp;
}

void writeVectorToField(fvMesh& mesh, const word& fld, const Eigen::VectorXd& y)
{
    auto makeNew = [&]()
    {
        return volVectorField
        (
            IOobject(fld, mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
            mesh,
            dimensionedVector(fld, dimsOf(fld, /*isVector*/true), vector::zero)
        );
    };

    volVectorField* vp = nullptr;
    {
        IOobject hdr(fld, "0", mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE);
        vp = hdr.typeHeaderOk<volVectorField>(true) ? new volVectorField(hdr, mesh)
                                                    : new volVectorField(makeNew());
    }
    auto& v = *vp;

    if (static_cast<label>(y.size()) != nEntriesWithBoundary(v))
        FatalErrorInFunction << "Size mismatch writing " << fld << exit(FatalError);

    label nCell = v.internalField().size();
    label k = 0;
    for (label i=0;i<nCell;++i)                        // internal
    {
        v[i].x() = y(k++); v[i].y() = y(k++); v[i].z() = y(k++);
    }

    auto& bf = v.boundaryFieldRef();                   // boundary (non-processor)
    forAll(bf, pi)
        if (!isA<processorFvPatch>(bf[pi].patch()))
            forAll(bf[pi], f)
            {
                bf[pi][f].x() = y(k++);
                bf[pi][f].y() = y(k++);
                bf[pi][f].z() = y(k++);
            }

    updateOnlyProcessorBC(v);                          // processor 동기화만
    v.write();
    delete vp;
}

} // namespace

// ------------------------------------------------------------

int main(int argc, char *argv[])
{
    argList::addBoolOption("noFunctionObjects", "do not execute functionObjects");

    #include "setRootCase.H"
    #include "createTime.H"

    // ---- 폴더 스캔으로 멀티/싱글 판정 & base mesh 생성 ----
    wordList regions; bool multiRegion=false;
    scanRegions(runTime, regions, multiRegion);

    const word baseRegion = multiRegion ? regions[0] : fvMesh::defaultRegion;

    fvMesh mesh
    (
        IOobject(baseRegion, runTime.timeName(), runTime, IOobject::MUST_READ)
    );

    if (!multiRegion) regions = wordList(1, mesh.name());

    Info << "regions: " << regions << nl << endl;

    // ---- fields ----
    const fileName listFile = "listField.dat";
    DynamicList<FieldSpec> all = readFieldList(listFile);

    // ---- 입력 매트릭스(Vinput, CurrentInput) ----
    const fileName baseDir = "constant";
    MatrixXd Vinput_global = readMatrixFlexible(baseDir/"Vinput.dat");     // (nSamples x nDim)
    MatrixXd CurrentInputM = readMatrixFlexible(baseDir/"CurrentInput.dat"); // (nDim x 1) 또는 (1 x nDim)
    VectorXd CurrentInput;
    if (CurrentInputM.cols() == 1) {
        CurrentInput = CurrentInputM.col(0);
    } else if (CurrentInputM.rows() == 1) {
        CurrentInput = CurrentInputM.row(0).transpose();
    } else {
        FatalErrorInFunction
            << "CurrentInput must be a vector: got "
            << CurrentInputM.rows() << "x" << CurrentInputM.cols() << exit(FatalError);
    }
    const int nDim = static_cast<int>(Vinput_global.cols());
    if (CurrentInput.size() != nDim)
        FatalErrorInFunction << "Dimension mismatch: CurrentInput("<<CurrentInput.size()
                             << ") vs Vinput cols("<<nDim<<")" << exit(FatalError);

    const int myRank = Pstream::myProcNo();

    for (const word& rName : regions)
    {
        fvMesh* rMeshPtr = (rName==mesh.name())
            ? &mesh
            : new fvMesh(IOobject(rName, runTime.timeName(), runTime, IOobject::MUST_READ));
        fvMesh& rMesh = *rMeshPtr;

        // region 필드 선택
        DynamicList<FieldSpec> myFields;
        forAll(all, i)
        {
            const auto& fs = all[i];
            if (fs.region == "-" || fs.region == rName) myFields.append(fs);
        }
        if (myFields.empty())
        {
            Info<< "  No fields assigned for region " << rName
                << " (check listField.dat)." << nl;
            if (rMeshPtr != &mesh) delete rMeshPtr;
            continue;
        }

        forAll(myFields, fi)
        {
            const word& fieldName = myFields[fi].name;
            const bool vecField   = isVector(myFields[fi].type);

            // --- POD 결과 읽기 (ASCII 헤더 형식)
            const fileName baseDir = "constant";
            const fileName fnUB   = baseDir/("UBasis_"+rName+"_"+fieldName+"_processor"+name(myRank)+".dat");
            const fileName fnA    = baseDir/("Acoef_"+rName+"_"+fieldName+"_processor"+name(myRank)+".dat");
            const fileName fnMean = baseDir/("MeanField_"+rName+"_"+fieldName+"_processor"+name(myRank)+".dat");

            MatrixXd UBasis = readMatrixFlexible(fnUB);      // (nX x nCase)
            MatrixXd Acoef  = readMatrixFlexible(fnA);       // (nCase x nCase): row = sample
            VectorXd mean   = readVectorFlexible(fnMean);    // (nX)

            const int nXcoord = static_cast<int>(UBasis.rows());
            const int nCase   = static_cast<int>(UBasis.cols());
            if (Acoef.rows() != nCase)
                FatalErrorInFunction << "Acoef rows("<<Acoef.rows()<<") != nCase("<<nCase<<") for field "<<fieldName
                                     << " region "<<rName << exit(FatalError);

            // --- nNearest = min(4, nCase, nSamplesAvail)
            const int nSamplesAvail = static_cast<int>(Vinput_global.rows());
            int nNearest = std::min(4, nCase);
            nNearest = std::min(nNearest, nSamplesAvail);
            if (nNearest < 1)
                FatalErrorInFunction << "No samples available for interpolation." << exit(FatalError);

            // --- 평가 벡터 [x, 1]
            MatrixXd EvaluationSample(nDim+1, 1);
            for (int iDim=0; iDim<nDim; ++iDim) EvaluationSample(iDim,0) = CurrentInput(iDim);
            EvaluationSample(nDim,0) = 1.0;

            // --- 최근접 nNearest 샘플 선별
            std::vector<double> distList(static_cast<size_t>(nSamplesAvail), 0.0);
            for (int iSample=0; iSample<nSamplesAvail; ++iSample)
            {
                double d=0.0;
                for (int j=0;j<nDim;++j)
                {
                    double diff = Vinput_global(iSample,j) - CurrentInput(j);
                    d += diff*diff;
                }
                distList[static_cast<size_t>(iSample)] = std::sqrt(d);
            }
            std::vector<size_t> idx(distList.size());
            std::iota(idx.begin(), idx.end(), 0);
            std::partial_sort(idx.begin(), idx.begin()+nNearest, idx.end(),
                              [&](size_t a, size_t b){ return distList[a] < distList[b]; });

            // --- 보간 시스템 구성
            MatrixXd CoordMat(nNearest, nDim+1); // [X | 1]
            MatrixXd FuncMat(nNearest, nCase);   // Acoef row 추출
            for (int i=0;i<nNearest;++i)
            {
                int id = static_cast<int>(idx[static_cast<size_t>(i)]);
                for (int j=0;j<nDim;++j) CoordMat(i,j) = Vinput_global(id,j);
                CoordMat(i,nDim) = 1.0;
                FuncMat.row(i)   = Acoef.row(id);
            }

            // --- 선형 회귀 보간 or 저차 보간
            MatrixXd InterpolatedACoef(nCase, 1);
            InterpolatedACoef.setZero();

            if (nNearest >= nDim + 1)
            {
                // 충분한 샘플: 최소제곱 LSQ
                Info << "LSQ interpolation: field " << fieldName << " region " << rName
                     << " (nNearest=" << nNearest << ", nDim=" << nDim << ")" << nl;
                MatrixXd coeff = CoordMat.colPivHouseholderQr().solve(FuncMat); // (nDim+1 x nCase)
                InterpolatedACoef = EvaluationSample.transpose() * coeff;       // (1 x nCase)
                InterpolatedACoef.transposeInPlace();                           // -> (nCase x 1)
            }
            else
            {
                // 샘플 부족: 차원 축소 보간 (SVD)
                Info << "Reduced-dim interpolation: field " << fieldName << " region " << rName
                     << " (nNearest=" << nNearest << ", nDim=" << nDim << ")" << nl;

                MatrixXd X(nNearest, nDim);
                for (int i=0;i<nNearest;++i)
                {
                    int id = static_cast<int>(idx[static_cast<size_t>(i)]);
                    for (int j=0;j<nDim;++j) X(i,j) = Vinput_global(id,j);
                }
                RowVectorXd mu = X.colwise().mean();
                MatrixXd Xc = X.rowwise() - mu;

                JacobiSVD<MatrixXd> svd(Xc, ComputeThinU | ComputeThinV);
                int rank = 0;
                {
                    const auto& S = svd.singularValues();
                    const double eps = 1e-12;
                    for (int k=0;k<S.size();++k) if (S(k) > eps) ++rank;
                }
                int r = std::min(std::min(rank, nDim), std::max(0, nNearest-1));

                if (r <= 0)
                {
                    // 상수 모델: FuncMat 행 평균
                    RowVectorXd rowMean = FuncMat.colwise().mean();
                    InterpolatedACoef = rowMean.transpose();
                }
                else
                {
                    MatrixXd Vr = svd.matrixV().leftCols(r); // (nDim x r)
                    MatrixXd Z  = Xc * Vr;                   // (nNearest x r)

                    MatrixXd CoordR(nNearest, r+1);
                    CoordR.leftCols(r) = Z;
                    CoordR.col(r).setOnes();

                    MatrixXd coeff_r = CoordR.colPivHouseholderQr().solve(FuncMat); // (r+1 x nCase)

                    RowVectorXd xEval(nDim);
                    for (int j=0;j<nDim;++j) xEval(j) = CurrentInput(j);
                    RowVectorXd zEval = (xEval - mu) * Vr;                           // (1 x r)

                    MatrixXd EvalR(r+1,1);
                    for (int j=0;j<r;++j) EvalR(j,0) = zEval(j);
                    EvalR(r,0) = 1.0;

                    MatrixXd pred = EvalR.transpose() * coeff_r; // (1 x nCase)
                    InterpolatedACoef = pred.transpose();        // (nCase x 1)
                }
            }

            // --- Reconstruction: UBasis * a + mean
            if (InterpolatedACoef.rows() != nCase)
                FatalErrorInFunction << "Interpolated a size mismatch: "
                                     << InterpolatedACoef.rows() << " vs nCase " << nCase << exit(FatalError);

            VectorXd y = UBasis * InterpolatedACoef; // (nX)
            if (mean.size() == nXcoord) y += mean;
            else Info << "Mean field missing or size mismatch; skipping mean add-back" << nl;

            // --- 출력 time
            runTime.setTime(10001.0, 10001);
            if (vecField) writeVectorToField(rMesh, fieldName, y);
            else          writeScalarToField(rMesh, fieldName, y);

            Info<< "  Reconstructed: region="<<rName<<", field="<<fieldName
                << " (nX="<<nXcoord<<", nCase="<<nCase<<", nNearest="<<nNearest<<")" << nl;
        }

        if (rMeshPtr != &mesh) delete rMeshPtr;
    }

    Info << "\n[baramPODreconstruct] Completed." << nl << endl;
    return 0;
}
