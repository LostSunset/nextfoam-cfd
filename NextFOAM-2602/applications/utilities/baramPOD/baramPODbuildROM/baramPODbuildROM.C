/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | baramPODbuildROM - folder-scan only
   \\    /   O peration     | listField.dat: <field> <type> <regionName|->
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"
#include "IOmanip.H"
#include "IFstream.H"
#include "OFstream.H"
#include "processorFvPatch.H"

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

// listField.dat
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

// 존재 체크
static inline bool scalarFieldExists(const fvMesh& mesh, const word& timeName, const word& fld)
{
    IOobject io(fld, timeName, mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE);
    return io.typeHeaderOk<volScalarField>(false);
}
static inline bool vectorFieldExists(const fvMesh& mesh, const word& timeName, const word& fld)
{
    IOobject io(fld, timeName, mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE);
    return io.typeHeaderOk<volVectorField>(false);
}

// MUST_READ 로 로드
template<class T>
auto readStrictField(const fvMesh& mesh, const word& timeName, const word& fld)
-> tmp<GeometricField<T, fvPatchField, volMesh>>
{
    IOobject io(fld, timeName, mesh, IOobject::MUST_READ, IOobject::NO_WRITE);
    if (!io.typeHeaderOk<GeometricField<T, fvPatchField, volMesh>>(false))
    {
        FatalErrorInFunction << "Missing field " << fld
                             << " at time " << timeName
                             << " for mesh " << mesh.name() << exit(FatalError);
    }
    return tmp<GeometricField<T, fvPatchField, volMesh>>(new GeometricField<T, fvPatchField, volMesh>(io, mesh));
}

VectorXd flattenField(const volScalarField& s)
{
    label n = s.internalField().size();
    const auto& bf = s.boundaryField();
    forAll(bf, pi)
        if (!isA<processorFvPatch>(bf[pi].patch()))
            n += bf[pi].size();

    VectorXd v(n);
    label k = 0;
    forAll(s, i) v(k++) = s[i];
    forAll(bf, pi)
        if (!isA<processorFvPatch>(bf[pi].patch()))
            forAll(bf[pi], f) v(k++) = bf[pi][f];
    return v;
}

VectorXd flattenField(const volVectorField& vv)
{
    label nCell = vv.internalField().size();
    label n = nCell;
    const auto& bf = vv.boundaryField();
    forAll(bf, pi)
        if (!isA<processorFvPatch>(bf[pi].patch()))
            n += bf[pi].size();

    VectorXd v(3*n);
    // internal
    label k = 0;
    for (label i=0;i<nCell;++i)
    {
        const vector& a = vv[i];
        v(k++) = a.x(); v(k++) = a.y(); v(k++) = a.z();
    }
    // boundary (non-processor)
    forAll(bf, pi)
        if (!isA<processorFvPatch>(bf[pi].patch()))
            forAll(bf[pi], f)
            {
                const vector& a = bf[pi][f];
                v(k++) = a.x(); v(k++) = a.y(); v(k++) = a.z();
            }
    return v;
}

// Eigen I/O
void writeEigen(const fileName& fn, const MatrixXd& M)
{
    OFstream os(fn);
    if (!os.good())
        FatalErrorInFunction << "Cannot write: " << fn << exit(FatalError);

    os  << M.rows() << " " << M.cols() << "\n";
    for (Eigen::Index i=0;i<M.rows();++i)
    {
        for (Eigen::Index j=0;j<M.cols();++j)
        {
            os  << M(i,j);
            if (j+1<M.cols()) os << " ";
        }
        os << "\n";
    }
}

void writeVector(const fileName& fn, const VectorXd& v)
{
    MatrixXd M(v.size(), 1);
    M.col(0) = v;
    writeEigen(fn, M);
}

// --- 폴더 스캔: constant/*/polyMesh/points
static void scanRegions(const Time& runTime, wordList& regions, bool& multiRegion)
{
    regions.clear();
    const fileName constDir = runTime.constant();

    fileNameList kids = readDir(constDir, fileName::DIRECTORY);
    forAll(kids, i)
    {
        const fileName sub = constDir/kids[i];
        const fileName pointsFile = sub/"polyMesh"/"points";
        if (isFile(pointsFile)) regions.append(kids[i]);
    }
    multiRegion = !regions.empty();
}

// 바꿔치기: 컨테이너 통째 reduce 대신, 원소별 스칼라 reduce
static void allreduceMatrixSum(Eigen::MatrixXd& M)
{
    // M.data()는 row-major가 아니어도 연속 메모리(Eigen default: column-major)
    // 원소 개수만큼 reduce를 한 번씩 수행 (간단, 호환성 ↑)
    for (Eigen::Index i = 0; i < M.size(); ++i)
    {
        Foam::scalar v = M.data()[i];
        Foam::reduce(v, Foam::sumOp<Foam::scalar>());  // 스칼라 합산 Allreduce
        M.data()[i] = v;
    }
}

} // namespace

// --------------------------------------------------------------

int main(int argc, char *argv[])
{
    argList::addBoolOption("noFunctionObjects", "do not execute functionObjects");
    timeSelector::addOptions(true);

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    if (timeDirs.empty())
    {
        FatalErrorInFunction << "No time directories found." << exit(FatalError);
    }

    // ---- 폴더 스캔으로 멀티/싱글 판정 & base mesh 생성 ----
    wordList regions; bool multiRegion=false;
    scanRegions(runTime, regions, multiRegion);

    const word baseRegion = multiRegion ? regions[0] : fvMesh::defaultRegion;

    fvMesh mesh
    (
        IOobject(baseRegion, runTime.timeName(), runTime, IOobject::MUST_READ)
    );

    // 최종 region 리스트 구성
    if (!multiRegion) regions = wordList(1, mesh.name());

    Info << "regions: " << regions << nl << endl;

    // ---- fields ----
    const fileName listFile = "listField.dat";
    DynamicList<FieldSpec> all = readFieldList(listFile);

    // 사용할 time (최신부터 역순)
    DynamicList<word> useTimes;
    forAll(timeDirs, ti) useTimes.append(timeDirs[ti].name());
    if (useTimes.empty())
        FatalErrorInFunction << "No usable time directories." << exit(FatalError);

    for (const word& rName : regions)
    {
        fvMesh* rMeshPtr = (rName==mesh.name())
            ? &mesh
            : new fvMesh(IOobject(rName, runTime.timeName(), runTime, IOobject::MUST_READ));
        fvMesh& rMesh = *rMeshPtr;

        // region용 필드 선택 (- 또는 정확 일치)
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
            if (rMeshPtr!=&mesh) delete rMeshPtr;
            continue;
        }

        // 각 필드 처리
        forAll(myFields, fi)
        {
            const word& fld = myFields[fi].name;
            const bool vec  = isVector(myFields[fi].type);

            DynamicList<VectorXd> snapshots;
            label nX = -1;

            forAll(useTimes, ti)
            {
                const word& tName = useTimes[ti];

                bool ok = vec ? vectorFieldExists(rMesh, tName, fld)
                              : scalarFieldExists(rMesh, tName, fld);
                if (!ok) continue;

                if (tName == "0") continue;  // 초기조건 타임스텝은 스킵
                if (tName == "10001") continue;

                VectorXd ref = vec ?
                    flattenField(readStrictField<vector>(rMesh, useTimes[1], fld)()) :
                    flattenField(readStrictField<scalar>(rMesh, useTimes[1], fld)());
                if (nX<0) nX = ref.size();

                VectorXd flat = vec ?
                    flattenField(readStrictField<vector>(rMesh, tName, fld)()) :
                    flattenField(readStrictField<scalar>(rMesh, tName, fld)());

                if (flat.size()!=nX)
                    FatalErrorInFunction << "Size mismatch for field " << fld
                                         << " at time " << tName << " in region " << rName
                                         << " (expected n=" << nX << ", got " << flat.size() << ")"
                                         << exit(FatalError);
                snapshots.append(std::move(flat));
            }

            if (snapshots.empty())
            {
                Info<< "  Skip (no snapshots): region="<<rName<<" field="<<fld<<nl;
                continue;
            }

            const label nCase = snapshots.size();
            MatrixXd Y(nX, nCase);
            for (int j=0;j<nCase;++j) Y.col(j) = snapshots[j];

            VectorXd mean = Y.rowwise().mean();
            MatrixXd Yc   = Y.colwise() - mean;

            MatrixXd C = (Yc.transpose()*Yc).eval();
            allreduceMatrixSum(C);
            SelfAdjointEigenSolver<MatrixXd> es(C);
            VectorXd L = es.eigenvalues();
            MatrixXd V = es.eigenvectors();

            // 내림차순 정렬
            {
                std::vector<std::pair<double,int>> order(L.size());
                for (int i=0;i<L.size();++i) order[i] = {L(i), i};
                std::sort(order.begin(), order.end(), [](auto&a, auto&b){return a.first>b.first;});
                VectorXd Ls(L.size()); MatrixXd Vs(V.rows(), V.cols());
                for (int k = 0; k < static_cast<int>(order.size()); ++k){ Ls(k)=order[k].first; Vs.col(k)=V.col(order[k].second); }
                L.swap(Ls); V.swap(Vs);
            }

            MatrixXd UB(nX, nCase); UB.setZero();
            for (int j=0;j<nCase;++j)
            {
                double lam=std::max(L(j),0.0);
                if(lam>SMALL) UB.col(j)=(Yc*V.col(j))/std::sqrt(lam);
            }

            MatrixXd Acoef(nCase,nCase);
            Acoef.setZero();
            for (int j=0; j<nCase; ++j)
            {
                const double sj = std::sqrt(std::max(L(j), 0.0));
                if (sj > SMALL) Acoef.col(j) = sj * V.col(j);
            }

            const int rank=Pstream::myProcNo();
            const fileName baseDir="constant";
            writeEigen(baseDir/("UBasis_"+rName+"_"+fld+"_processor"+Foam::name(rank)+".dat"),UB);
            writeEigen(baseDir/("Acoef_"+rName+"_"+fld+"_processor"+Foam::name(rank)+".dat"),Acoef);
            writeEigen(baseDir/("EigenValues_"+rName+"_"+fld+"_processor"+Foam::name(rank)+".dat"),L);
            writeVector(baseDir/("MeanField_"+rName+"_"+fld+"_processor"+Foam::name(rank)+".dat"),mean);

            Info<<"  POD done: region="<<rName<<" field="<<fld<<nl;
        }

        if (rMeshPtr!=&mesh) delete rMeshPtr;
    }

    Info<<"\n[baramPODbuildROM] Completed."<<nl;
    return 0;
}
