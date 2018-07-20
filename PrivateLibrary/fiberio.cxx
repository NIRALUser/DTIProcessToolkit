#include <string>
#include <cmath>
#include <memory>

#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectWriter.h>
#include <vtkVersion.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkCell.h>
#include <vtkFloatArray.h>

#include "fiberio.h"

// hide function to this compilation unit
namespace
{
inline double SQ2(double x)
{
  return x * x;
}

};

void writeFiberFile(const std::string & filename, GroupType::Pointer fibergroup, bool saveProperties, std::string encoding, std::string scalarPropertyName)
{
  // Make sure origins are updated
  fibergroup->ComputeObjectToWorldTransform();

  // ITK Spatial Object
  if( filename.rfind(".fib") != std::string::npos )
    {
    typedef itk::SpatialObjectWriter<3> WriterType;
    WriterType::Pointer writer  = WriterType::New();
    writer->SetInput(fibergroup);
    writer->SetFileName(filename);
    writer->Update();
    }
  // VTK Poly Data
  else if( filename.rfind(".vt") != std::string::npos )
    {
    // Build VTK data structure
    vtkSmartPointer<vtkPolyData>   polydata(vtkPolyData::New() );
    vtkSmartPointer<vtkFloatArray> tensorsdata(vtkFloatArray::New() );
    vtkSmartPointer<vtkIdList>     ids(vtkIdList::New() );
    vtkSmartPointer<vtkPoints>     pts(vtkPoints::New() );

    tensorsdata->SetNumberOfComponents(9);
    polydata->SetPoints(pts);

    ids->SetNumberOfIds(0);
    pts->SetNumberOfPoints(0);
    polydata->Allocate();
    vtkFloatArray *scalarFA = vtkFloatArray::New();
    vtkFloatArray *scalarMD = vtkFloatArray::New();
    vtkFloatArray *scalarRD = vtkFloatArray::New();
    vtkFloatArray *scalarAD = vtkFloatArray::New();
    vtkFloatArray *scalarValue = vtkFloatArray::New();
    scalarFA->SetNumberOfComponents(1);
    scalarFA->SetName("FA");
    scalarMD->SetNumberOfComponents(1);
    scalarMD->SetName("MD");
    scalarAD->SetNumberOfComponents(1);
    scalarAD->SetName("AD");
    scalarRD->SetNumberOfComponents(1);
    scalarRD->SetName("RD");
    scalarValue->SetNumberOfComponents(1);
    if (scalarPropertyName != "")
      {
	scalarValue->SetName(scalarPropertyName.c_str());
      }

  

    std::auto_ptr<ChildrenListType> children(fibergroup->GetChildren(0) );
    typedef ChildrenListType::const_iterator IteratorType;
    for( IteratorType it = children->begin(); it != children->end(); it++ )
      {
      itk::SpatialObject<3>*        tmp = (*it).GetPointer();
      itk::DTITubeSpatialObject<3>* tube = dynamic_cast<itk::DTITubeSpatialObject<3> *>(tmp);
      unsigned int                  nPointsOnFiber = tube->GetNumberOfPoints();
      vtkIdType                     currentId = ids->GetNumberOfIds();
      for( unsigned int k = 0; k < nPointsOnFiber; k++ )
        {
        itk::Point<double, 3>  v(tube->GetPoint(k)->GetPosition() );
        itk::Vector<double, 3> spacing(tube->GetSpacing() );
        itk::Vector<double, 3> origin(tube->GetObjectToWorldTransform()->GetOffset() );

        // convert origin from LPS -> RAS
        origin[0] = -origin[0];
        origin[1] = -origin[1];

        vtkIdType id;
        // Need to multiply v by spacing and origin
        // Also negate the first to convert from LPS -> RAS
        // for slicer 3
        id = pts->InsertNextPoint(-v[0] * spacing[0] + origin[0],
                                  -v[1] * spacing[1] + origin[1],
                                  v[2] * spacing[2] + origin[2]);

        ids->InsertNextId(id);

        itk::DTITubeSpatialObjectPoint<3>* sopt = dynamic_cast<itk::DTITubeSpatialObjectPoint<3> *>(tube->GetPoint(k) );
        float                              vtktensor[9];
        vtktensor[0] = sopt->GetTensorMatrix()[0];
        vtktensor[1] = sopt->GetTensorMatrix()[1];
        vtktensor[2] = sopt->GetTensorMatrix()[2];
        vtktensor[3] = sopt->GetTensorMatrix()[1];
        vtktensor[4] = sopt->GetTensorMatrix()[3];
        vtktensor[5] = sopt->GetTensorMatrix()[4];
        vtktensor[6] = sopt->GetTensorMatrix()[2];
        vtktensor[7] = sopt->GetTensorMatrix()[4];
        vtktensor[8] = sopt->GetTensorMatrix()[5];

	//std::cout << "MD = " << sopt->GetField("md") << ", " ;

        scalarFA->InsertNextValue(sopt->GetField("fa"));
        scalarMD->InsertNextValue(sopt->GetField("md"));
        scalarAD->InsertNextValue(sopt->GetField("ad"));
        scalarRD->InsertNextValue(sopt->GetField("rd"));


#if VTK_MAJOR_VERSION > 7 || ( VTK_MAJOR_VERSION == 7 && VTK_MINOR_VERSION >=1 )
        tensorsdata->InsertNextTypedTuple(vtktensor);
#else
        tensorsdata->InsertNextTupleValue(vtktensor);
#endif

	if (scalarPropertyName != "")
	  {
	    scalarValue->InsertNextValue(sopt->GetField(scalarPropertyName.c_str()));
	  }

        }
      polydata->InsertNextCell(VTK_POLY_LINE, nPointsOnFiber, ids->GetPointer(currentId) );
      }

    polydata->GetPointData()->SetTensors(tensorsdata);
    if (saveProperties)
      {
	polydata->GetPointData()->AddArray(scalarFA);
	polydata->GetPointData()->AddArray(scalarMD);
	polydata->GetPointData()->AddArray(scalarAD);
	polydata->GetPointData()->AddArray(scalarRD);
      }
    if (scalarPropertyName != "")
      {
	polydata->GetPointData()->AddArray(scalarValue);
      }

    // Legacy
    if( filename.rfind(".vtk") != std::string::npos )
      {
      vtkSmartPointer<vtkPolyDataWriter> fiberwriter = vtkPolyDataWriter::New();
      fiberwriter->SetFileTypeToBinary();
      fiberwriter->SetFileName(filename.c_str() );
#if (VTK_MAJOR_VERSION < 6)
      fiberwriter->SetInput(polydata);
#else
      fiberwriter->SetInputData(polydata);
#endif
      if( encoding == "binary" )
      {
        fiberwriter->SetFileTypeToBinary();
      }
      else
      {
        fiberwriter->SetFileTypeToASCII() ;
      }
      fiberwriter->Update();
      }
    // XML
    else if( filename.rfind(".vtp") != std::string::npos )
      {
      vtkSmartPointer<vtkXMLPolyDataWriter> fiberwriter = vtkXMLPolyDataWriter::New();
      fiberwriter->SetFileName(filename.c_str() );
#if (VTK_MAJOR_VERSION < 6)
      fiberwriter->SetInput(polydata);
#else
      fiberwriter->SetInputData(polydata);
#endif
      if( encoding == "binary" )
      {
        fiberwriter->SetDataModeToBinary() ;
      }
      else if( encoding == "appended" )
      {
        fiberwriter->SetDataModeToAppended() ;
      }
      else
      {
        fiberwriter->SetDataModeToAscii() ;
      }
      fiberwriter->Update();
      }
    else
      {
      throw itk::ExceptionObject("Unknown file format for fibers");
      }
    }
  else
    {
    throw itk::ExceptionObject("Unknown file format for fibers");
    }
}

GroupType::Pointer readFiberFile(const std::string & filename)
{

  // ITK Spatial Object
  if( filename.rfind(".fib") != std::string::npos )
    {
    typedef itk::SpatialObjectReader<3, unsigned char> SpatialObjectReaderType;

    // Reading spatial object
    SpatialObjectReaderType::Pointer soreader = SpatialObjectReaderType::New();

    soreader->SetFileName(filename);
    soreader->Update();

    return soreader->GetGroup();
    }
  // VTK Poly Data
  else if( filename.rfind(".vt") != std::string::npos )
    {
    // Build up the principal data structure for fiber tracts
    GroupType::Pointer fibergroup = GroupType::New();

    vtkSmartPointer<vtkPolyData> fibdata(ITK_NULLPTR);

    // Legacy
    if( filename.rfind(".vtk") != std::string::npos )
      {
      vtkSmartPointer<vtkPolyDataReader> reader(vtkPolyDataReader::New() );
      reader->SetFileName(filename.c_str() );
      reader->Update();
      fibdata = reader->GetOutput();

      }
    else if( filename.rfind(".vtp") != std::string::npos )
      {
      vtkSmartPointer<vtkXMLPolyDataReader> reader(vtkXMLPolyDataReader::New() );
      reader->SetFileName(filename.c_str() );
      reader->Update();
      fibdata = reader->GetOutput();
      }
    else
      {
      throw itk::ExceptionObject("Unknown file format for fibers");
      }

    typedef  itk::SymmetricSecondRankTensor<double, 3> ITKTensorType;
    typedef  ITKTensorType::EigenValuesArrayType       LambdaArrayType;

    // Iterate over VTK data
    const int nfib = fibdata->GetNumberOfCells();
    int       pindex = -1;
    for( int i = 0; i < nfib; ++i )
      {
      itk::DTITubeSpatialObject<3>::Pointer dtiTube = itk::DTITubeSpatialObject<3>::New();
      vtkSmartPointer<vtkCell>              fib = fibdata->GetCell(i);

      vtkSmartPointer<vtkPoints> points = fib->GetPoints();

      typedef itk::DTITubeSpatialObjectPoint<3> DTIPointType;
      std::vector<DTIPointType> pointsToAdd;

      vtkSmartPointer<vtkDataArray> fibtensordata = fibdata->GetPointData()->GetTensors();
      for( int j = 0; j < points->GetNumberOfPoints(); ++j )
	{
	  ++pindex;

	  double * coordinates = points->GetPoint(j);
	  DTIPointType          pt;
	  // Convert from RAS to LPS for vtk
	  pt.SetPosition(-coordinates[0], -coordinates[1], coordinates[2]);
	  pt.SetRadius(0.5);
	  pt.SetColor(0.0, 1.0, 0.0);
	  double * vtktensor;
	  float                 floattensor[6];
	  ITKTensorType         itktensor;

	  if (fibtensordata) {
	    vtktensor = fibtensordata->GetTuple9(pindex);
	    floattensor[0] = itktensor[0] = vtktensor[0];
	    floattensor[1] = itktensor[1] = vtktensor[1];
	    floattensor[2] = itktensor[2] = vtktensor[2];
	    floattensor[3] = itktensor[3] = vtktensor[4];
	    floattensor[4] = itktensor[4] = vtktensor[5];
	    floattensor[5] = itktensor[5] = vtktensor[8];
	  } else {
	    vtktensor = new double[9];
	    floattensor[0] = itktensor[0] = vtktensor[0] = 1.0;
	    floattensor[1] = itktensor[1] = vtktensor[1] = 0.0;
	    floattensor[2] = itktensor[2] = vtktensor[2] = 0.0;
	    floattensor[3] = itktensor[3] = vtktensor[4] = 1.0;
	    floattensor[4] = itktensor[4] = vtktensor[5] = 0.0;
	    floattensor[5] = itktensor[5] = vtktensor[8] = 1.0;
	  }

	  pt.SetTensorMatrix(floattensor);

	  LambdaArrayType lambdas;

	  // Need to do do eigenanalysis of the tensor
	  itktensor.ComputeEigenValues(lambdas);

	  // FIXME: We should not be repeating this code here.  The code
	  // for all these computations should be re-factored into a
	  // common library.

	  float md = (lambdas[0] + lambdas[1] + lambdas[2]) / 3;
	  float fa = sqrt(1.5) * sqrt( (lambdas[0] - md) * (lambdas[0] - md)
				       + (lambdas[1] - md) * (lambdas[1] - md)
				       + (lambdas[2] - md) * (lambdas[2] - md) )
	    / sqrt(lambdas[0] * lambdas[0] + lambdas[1] * lambdas[1] + lambdas[2] * lambdas[2]);

	  float logavg = (log(lambdas[0]) + log(lambdas[1]) + log(lambdas[2]) ) / 3;

	  float ga =  sqrt( SQ2(log(lambdas[0]) - logavg) \
			    + SQ2(log(lambdas[1]) - logavg)	\
			    + SQ2(log(lambdas[2]) - logavg) );

	  float rd = (lambdas[1] + lambdas[0])/2;

	  pt.AddField("fa", fa);
	  pt.AddField("ga", ga);
	  pt.AddField("md", md);
	  pt.AddField("l1", lambdas[2]);
	  pt.AddField("l2", lambdas[1]);
	  pt.AddField("l3", lambdas[0]);
	  pt.AddField("rd", rd);

	  pointsToAdd.push_back(pt);
        }

      dtiTube->SetPoints(pointsToAdd);
      fibergroup->AddSpatialObject(dtiTube);
      }
    return fibergroup;
    } // end process .vtk .vtp
  else
    {
    throw itk::ExceptionObject("Unknown fiber file");
    }
}
