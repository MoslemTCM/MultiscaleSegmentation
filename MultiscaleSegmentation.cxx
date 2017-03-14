/*========================================================================= */

// Équipe Télédetection pour les catastrophes majeures (TCM).

// Programme : Application d'une segmentation multiechelle basée sur la transformée en ondelettes et watershed.

// Auteur : Moslem Ouled Sghaier

// Version : 0

/*========================================================================= */
#include <stdlib.h>

#include "itkBinaryImageToLabelMapFilter.h"
#include "itkBinaryThinningImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkChangeInformationImageFilter.h"
#include "itkLabelContourImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLabelImageToLabelMapFilter.h" 
#include "itkLabelMap.h"
#include "itkLabelMapToBinaryImageFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelObject.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkScalarToRGBPixelFunctor.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkWatershedImageFilter.h"


#include "otbAttributesMapLabelObject.h"
#include "otbBandMathImageFilter.h"
#include "otbImage.h"
#include "otbImageFileReader.h"
#include "otbImageFileWriter.h"
#include "otbLabelImageToLabelMapWithAdjacencyFilter.h"
#include "otbScalarImageToTexturesFilter.h"
#include "otbScalarImageToAdvancedTexturesFilter.h"
#include "otbShapeAttributesLabelMapFilter.h"
#include "otbStatisticsAttributesLabelMapFilter.h"
#include "otbStreamingImageFileWriter.h"
#include "otbWatershedSegmentationFilter.h"
#include "otbWaveletFilterBank.h"
#include "otbWaveletGenerator.h"
#include "otbWaveletOperator.h"
#include "otbWaveletTransform.h"

#include "otbWrapperApplication.h" 
#include "otbWrapperApplicationRegistry.h"
#include "otbWrapperApplicationFactory.h"
#include "otbWrapperTags.h"

//Utils
#include "itksys/SystemTools.hxx"
#include "itkListSample.h"

// Elevation handler
#include "otbWrapperElevationParametersHandler.h"


namespace otb
{
namespace Wrapper
{

class MultiscaleSegmentation : public Application
{

public:
	typedef MultiscaleSegmentation Self;
    typedef Application              Superclass;
    typedef itk::SmartPointer<Self>       Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

/** Standard macro */
    itkNewMacro(Self);
    itkTypeMacro(MultiscaleSegmentation, otb::Application);




// ********************************************************************************************************************
// Tuning parameters
// ********************************************************************************************************************

static const int labelObjectMinimumPixelSize = 15;
static const int textureRadius = 5;


// ********************************************************************************************************************


static const int Dimension = 2;
typedef float WaveletPixelType; // To our knowledge, Wavelet works only with "float".
typedef itk::Image<WaveletPixelType, 2>	WaveletImageType;


// ********************************************************************************************************************
// Prototypes

//int waveletTransformProcess(int argc, char* argv[]);
//int segmentationProcess(int argc, char * argv[]);

/* ********************************************************************************************************************

Fonction d'entrée

******************************************************************************************************************** */

private:


void DoInit()
{

SetName("MultiscaleSegmentation"); // Nécessaire
SetDocName("MultiscaleSegmentation");
SetDocLongDescription("Un simple module pour la segmentation de l'image basée sur la transformée en ondelettes et watershed");
SetDocLimitations("Les autres paramètres seront ajoutés plus tard");
SetDocAuthors("Moslem Ouled Sghaier");


AddParameter(ParameterType_InputImage,"in", "Input Image");
SetParameterDescription("in", "The input image");
AddParameter(ParameterType_OutputImage,"out", "Output Image");
SetParameterDescription("out","The output image");

}

void DoUpdateParameters()
{
	// Nothing to do here : all parameters are independent
}


void DoExecute()

{

	std::cout << "Multiechelle : HAAR wavelet-watershed transform\n";
	std::cout << " Les arguments sont donnés par moi \n";
	std::cout << "ARG1 : Input file\n";
	std::cout << "ARG2 : Output file\n";
	std::cout << "ARG3 : Watershed.Level\n";
	std::cout << "ARG4 : Watershed.Threshold\n";
	std::cout << "ARG5 : Testure fusion threshold\n";
	
	char *inputFile = "";
	char *outputFile = "";
	char *watershedLevel ="0.30";
	char *watershedThreshold ="0.001";
	char *textureFusionThreshold ="0.0012";

	float multiscaleContourRatioOfMaximumIntensity = 0.30f;
    char *fusionIterationCountString = "5";

	// *** Calcul la transformée en ondelettes ***
	const int waveletArgCount = 5;
	char * waveletArg[waveletArgCount] =
	{
		"",
		inputFile,
		"",
		"2",
		"2"
	};


	 waveletTransformProcess(waveletArgCount, waveletArg);
	
	// *** Calcul de la segmentation ***
	const int segmentationArgCount = 7;

	char * segmentation0Arg[segmentationArgCount] =
	{
		"",
		"",
		watershedLevel,
		watershedThreshold,
		"2", // Niveau de la décomposition en ondelette
		textureFusionThreshold,
		fusionIterationCountString
	};

	segmentationProcess(segmentationArgCount, segmentation0Arg);
   
	char * segmentation1Arg[segmentationArgCount] =
	{
		"",
		"",
		watershedLevel,
		watershedThreshold,
		"1", // Niveau de la décomposition en ondelette
		textureFusionThreshold,
		fusionIterationCountString
	};

	segmentationProcess(segmentationArgCount, segmentation1Arg);


	char * segmentation2Arg[segmentationArgCount] =
	{
		"",
		outputFile,
		watershedLevel,
		watershedThreshold,
		"0", // Niveau de la décomposition en ondelette
		textureFusionThreshold,
		fusionIterationCountString
	};

	segmentationProcess(segmentationArgCount, segmentation2Arg);
	
}


/* ********************************************************************************************************************

Calcul la transformée en ondelettes de l'entrée fournie.

******************************************************************************************************************** */
int waveletTransformProcess(int argc, char* argv[])
{
	float multiscaleContourRatioOfMaximumIntensity = 0.30f;

	std::cout << "* HAAR wavelet transform processing *\n";
	std::cout << "ARG1: " << " File input : " << argv[1] << "\n";
	std::cout << "ARG2: " << " File output - NOT USED\n";
	std::cout << "ARG3: " << " wavelet.level : " << atoi(argv[3]) << "\n";
	std::cout << "ARG4: " << " wavelet.decimator : " << atoi(argv[4]) << "\n";
	std::cout << "\n" << argc - 1 << " argument provided\n";

	if(argc != 5)
	{
		std::cerr << "Wrong number of parameters\n";
		return EXIT_FAILURE;
	}

	const char *       inputFileName = argv[1];
	const char *       outputFileName = argv[2];
	const unsigned int requestedLevel = atoi(argv[3]);
	const unsigned int decimFactor = atoi(argv[4]);
	
	unsigned int output_startup_index = 0;
	if(argc >= 6 )
	{
		output_startup_index = atoi(argv[5]);
	}

	if (argc == 7)
	{
		unsigned int  NbOfThreads = atoi(argv[6]);
		itk::MultiThreader::SetGlobalDefaultNumberOfThreads(NbOfThreads);
	}

	std::cout << argc << " argument processed\n";


	typedef otb::Image<WaveletPixelType, Dimension> ImageType;

	ImageType::Pointer reader = GetParameterFloatImage("in");

	SetParameterOutputImage<ImageType>("out",reader);

	char output_name[50];  

	typedef otb::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();

	// Choix de l'ondelette
	const otb::Wavelet::Wavelet wvltID = otb::Wavelet::Wavelet::HAAR;

	// *** Applique la transformée ***
	typedef otb::WaveletOperator<wvltID, otb::Wavelet::FORWARD, WaveletPixelType, Dimension> WaveletOperator;
	typedef otb::WaveletFilterBank<ImageType, ImageType, WaveletOperator, otb::Wavelet::FORWARD>  ForwardFilterBank;
	typedef otb::WaveletTransform<ImageType, ImageType, ForwardFilterBank, otb::Wavelet::FORWARD> FilterType;


	FilterType::Pointer filter = FilterType::New();
	filter->SetInput(reader);

	for(unsigned int currentLevel = 1; currentLevel <= requestedLevel; ++currentLevel)
	{
		filter->SetNumberOfDecompositions(currentLevel);
		filter->SetSubsampleImageFactor(decimFactor);

		try
		{
			filter->Update();
		}
		catch (itk::ExceptionObject& err)
		{
			std::cout << "ExceptionObject caught !" << std::endl;
			std::cout << err << std::endl;
			return EXIT_FAILURE;
		}


		// *** Envoie la décomposition dans un fichier ***

		sprintf(output_name, "WaveletLevel%d.tif", currentLevel);

		writer->SetFileName(output_name);
		writer->SetInput(filter->GetOutput()->Front());

		try
		{
			writer->Update();
		}
		catch (itk::ExceptionObject& err)
		{
			std::cout << "ExceptionObject caught !" << std::endl;
			std::cout << err << std::endl;
			return EXIT_FAILURE;
		}

		std::cout<< output_name << " generated\n";

		// *** Envoie les différentes décomposition dans des fichiers ***

		int count = output_startup_index;

		for(auto image_iterator = filter->GetOutput()->Begin(); 
			image_iterator != filter->GetOutput()->End(); 
			++image_iterator)
		{
			sprintf(output_name, "WaveletMultiscaleOutput%d.tif",count);

			writer->SetFileName(output_name);
			writer->SetInput(image_iterator.Get());

			try
			{
				writer->Update();
			}
			catch (itk::ExceptionObject& err)
			{
				std::cout << "ExceptionObject caught !" << std::endl;
				std::cout << err << std::endl;
				return EXIT_FAILURE;
			}
		
			std::cout << output_name << " produced\n";

			++count;
		}

	}

	// *** Récupère l'image source ***
	writer->SetFileName("WaveletLevel0.tif");
	writer->SetInput(reader);

	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject& err)
	{
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "\n\nProcess completed.\n";


	return EXIT_SUCCESS;

}

/* ********************************************************************************************************************

Crée une image utilisée dans le redimensionnement d'image.

******************************************************************************************************************** */
void CreateImage(WaveletImageType::Pointer image)
{
	itk::Index<2> start; start.Fill(0);
	itk::Size<2> size; size.Fill(100);
	WaveletImageType::RegionType region(start, size);
	image->SetRegions(region);
	image->Allocate();
	image->FillBuffer(0);
 
	for(unsigned int r = 40; r < 60; r++)
	{
		for(unsigned int c = 40; c < 60; c++)
		{
			WaveletImageType::IndexType pixelIndex;
			pixelIndex[0] = r;
			pixelIndex[1] = c;
			image->SetPixel(pixelIndex, 255);
		}
	}
}

/* ********************************************************************************************************************

Imprime la liste des attributs d'un LabelObjet.

******************************************************************************************************************** */

static void PrintLabelObjectAttributes(otb::AttributesMapLabelObject<unsigned short, Dimension,	float>::Pointer thisLabelObject)
{

	if(thisLabelObject == nullptr)
	{
		return;
	}


	try
	{
		std::cout << "Number of attributes :" << thisLabelObject->GetNumberOfAttributes() << "\n";
	
		if(thisLabelObject->GetNumberOfAttributes() != 0)
		{
			typedef std::vector<std::string> attributeVectorType;
			attributeVectorType availableAttributeVector = thisLabelObject->GetAvailableAttributes();

			for(std::vector<std::string>::iterator attributeIterator = availableAttributeVector.begin();
				attributeIterator != availableAttributeVector.end();
				++attributeIterator)
			{
				std::cout << "Available attributes: " << *attributeIterator << "\n";;
			}
		}
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}

}

/* ********************************************************************************************************************

Segmente le niveau d'ondelette précisé.

******************************************************************************************************************** */

int segmentationProcess(int argc, char * argv[])
{
	float multiscaleContourRatioOfMaximumIntensity = 0.30f;

	std::cout << "*** Segmentation processing ***\n";
	std::cout << "ARG1: " << " Output file : " << argv[1] << "\n";
	std::cout << "ARG2: " << " watershed.level : " << atof(argv[2])<< "\n";
	std::cout << "ARG3: " << " watershed.threshold : " << atof(argv[3]) << "\n";
	std::cout << "ARG4: " << " decimation level : " << atoi(argv[4]) << "\n";
	std::cout << "ARG5: " << " Texture fusion threshold : " << atof(argv[5]) << "\n";
	std::cout << "ARG6: " << " Fusion iteration : " << atoi(argv[6]) << "\n";

	std::cout << "\n";

	if (argc != 7) 
	{
		std::cerr << "\nWrong number of parameters\n";
		return EXIT_FAILURE;
	}
	
	// ****************************************************************************************************************
	
	typedef unsigned char CharPixelType; // IO
    const unsigned int Dimension1 = 2;
    typedef otb::Image<CharPixelType, Dimension1> CharImageType;


	typedef itk::RGBPixel<unsigned char>				RGBPixelType;
	
	typedef otb::Image<RGBPixelType, 2>					RGBImageType;
	typedef itk::Image<unsigned long, 2>				WatershedLabeledImageType;

	typedef unsigned short                              LabelType;
	typedef otb::Image<LabelType, Dimension>            LabeledImageType;
	
	typedef otb::AttributesMapLabelObject<LabelType, Dimension,	float>		LabelObjectType;
	typedef otb::LabelMapWithAdjacency<LabelObjectType>						LabelMapType;
	typedef itk::LabelImageToLabelMapFilter<LabeledImageType, LabelMapType>	LabelMapFilterType;


	typedef otb::LabelImageToLabelMapWithAdjacencyFilter<WaveletImageType, LabelMapType> ImageToLabelMapFilterType;
	typedef itk::LabelMapToLabelImageFilter<LabelMapType, WaveletImageType> LabelMapToLabelImageFilterType;	

	typedef otb::StatisticsAttributesLabelMapFilter<LabelMapType,  WaveletImageType> StatisticsLabelMapFilterType;
		
	typedef otb::ImageFileReader<WaveletImageType> WatershedFileReaderType;	

	typedef otb::ImageFileWriter<WaveletImageType> postWriterType;
	typedef otb::ImageFileWriter<RGBImageType> postColorWriterType;

	typedef itk::RescaleIntensityImageFilter<WaveletImageType, WaveletImageType> RescaleIntensityImageFilterType;
	typedef itk::RescaleIntensityImageFilter<WaveletImageType, CharImageType> RescaleIntensityImageFilterType1;

	typedef itk::GradientMagnitudeImageFilter<WaveletImageType, WaveletImageType> GradientMagnitudeFilterType;


	typedef itk::LabelContourImageFilter<WaveletImageType, WaveletImageType> LabelImageContourFilterType;

	typedef itk::BinaryThinningImageFilter<WaveletImageType, WaveletImageType> BinaryThinningImageFilterType;

	typedef otb::BandMathImageFilter<WaveletImageType> BandMathFilterType;

	// ****************************************************************************************************************

	char filenameBuffer[200];

	int decimationLevel = atoi(argv[4]);
	
	unsigned int mergeIterationCount = atoi(argv[6]);

	float mergingThreshold = atof(argv[5]);

	char inputFileName[50];

	sprintf(inputFileName, "WaveletLevel%d.tif", decimationLevel);


	char labelImageFileName[50];
		
	switch(decimationLevel)
	{
		case 0:
			sprintf(labelImageFileName, "Level1LabelImage.tif");
			break;
		case 1:
			sprintf(labelImageFileName, "Level2LabelImage.tif");
			break;

		default:
		case 2:
			labelImageFileName[0] = 0;
			break;

	}

	char outputFileName[50];

	switch(decimationLevel)
	{
		case 0:
			sprintf(outputFileName, "Level0LabelImage.tif");
			break;
		case 1:
			sprintf(outputFileName, "Level1LabelImage.tif");
			break;
		default:
		case 2:
			sprintf(outputFileName, "Level2LabelImage.tif");
			break;

	}

	
	LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
	postWriterType::Pointer outputWriter = postWriterType::New();
	postColorWriterType::Pointer postColorWatershedWriter = postColorWriterType::New();
	
	
	// ****************************************************************************************************************
	// Étire la segmentation du niveau d'ondelettes précédent au niveau de l'ondelette courant
	// puis extrait les contours pour l'applicaton multiéchelle.
	

	LabelImageContourFilterType::Pointer edgeContourFilter = LabelImageContourFilterType::New();

	BinaryThinningImageFilterType::Pointer edgeThinnerFilter = BinaryThinningImageFilterType::New();

	typedef itk::IdentityTransform<double, 2> TransformType;
	typedef itk::ResampleImageFilter<WaveletImageType, WaveletImageType> ResampleImageFilterType;
	ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();

	if(decimationLevel < 2)
	{
		std::cout << "Stretching " << labelImageFileName << " .. \n";

		WatershedFileReaderType::Pointer StretchingReader = WatershedFileReaderType::New();

		StretchingReader->SetFileName(labelImageFileName);
 
		try
		{
			StretchingReader->Update();
		}
		catch (itk::ExceptionObject& e)
		{
			std::cerr << e << std::endl;
		}

		WaveletImageType::SizeType inputSize = StretchingReader->GetOutput()->GetLargestPossibleRegion().GetSize();

		std::cout << "    Input size: " << inputSize << std::endl;
			
		// Redemensionne
		WaveletImageType::SizeType outputSize;
		outputSize[0] = inputSize[0] * 2;
		outputSize[1] = inputSize[1] * 2;

		WaveletImageType::SpacingType outputSpacing;

		outputSpacing[0] = StretchingReader->GetOutput()->GetSpacing()[0] * (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
		outputSpacing[1] = StretchingReader->GetOutput()->GetSpacing()[1] * (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));
 

		resample->SetInput(StretchingReader->GetOutput());
		resample->SetSize(outputSize);
		resample->SetOutputSpacing(outputSpacing);
		resample->SetTransform(TransformType::New());
		resample->UpdateLargestPossibleRegion();
 
		WaveletImageType::Pointer output = resample->GetOutput();
 
		std::cout << "    Output size: " << output->GetLargestPossibleRegion().GetSize() << std::endl;

		sprintf(filenameBuffer, "Level%dLabelImageStretched.tif", decimationLevel + 1);

		outputWriter->SetFileName(filenameBuffer);
		outputWriter->SetInput(output);
		
		try
		{
			outputWriter->Update();
		}
		catch (itk::ExceptionObject& e)
		{
			std::cerr << e << std::endl;
		}

		// *** Effectue l'extraction du contour de l'image segmenté

		
		edgeContourFilter->SetInput(output);
		edgeContourFilter->SetFullyConnected(false); // Produce smaller contour

		try
		{
			edgeContourFilter->Update();
		}
		catch (itk::ExceptionObject& e)
		{
			std::cerr << e << std::endl;
		}

		edgeThinnerFilter->SetInput(edgeContourFilter->GetOutput());


		// *** Imprime les contours extrait sur l'image segmenté pour déboggage ***

		WatershedFileReaderType::Pointer contourPrintReader = WatershedFileReaderType::New();
		contourPrintReader->SetFileName(filenameBuffer);
		
		BandMathFilterType::Pointer contourOnStretched = BandMathFilterType::New();

		contourOnStretched->SetNthInput(0, contourPrintReader->GetOutput());

		contourOnStretched->SetNthInput(1, edgeThinnerFilter->GetOutput());
	

		contourOnStretched->SetExpression("if(b2 > 0, 500, b1)");

		sprintf(filenameBuffer,"Level%dLabelImageStretchedWithBoundary.tif", decimationLevel+1);
	 	outputWriter->SetFileName(filenameBuffer);
		outputWriter->SetInput(contourOnStretched->GetOutput());
		try
		{
			outputWriter->Update();
		}
		catch (itk::ExceptionObject& e)
		{
			std::cerr << e << std::endl;
		}

	}



	// ********************************************************************************************
	// *** Applique un gradient pré-watershed ***

	WatershedFileReaderType::Pointer inputReader = WatershedFileReaderType::New();
	inputReader->SetFileName(inputFileName);

	
	std::cout << "Computing gradient on " << inputFileName << " ..\n";


	GradientMagnitudeFilterType::Pointer inputGradient = GradientMagnitudeFilterType::New();
	inputGradient->SetUseImageSpacingOff();
	inputGradient->SetInput(inputReader->GetOutput());

	try
	{
		inputGradient->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}

	// ********************************************************************************************
	// *** Ajoute au pré-watershed le contour du niveau supérieur si nécessaire. ***

	BandMathFilterType::Pointer segmentationOnSource = BandMathFilterType::New();

	if(decimationLevel < 2)
	{
		
		segmentationOnSource->SetNthInput(0, inputGradient->GetOutput());
		
 
		typedef itk::ChangeInformationImageFilter< WaveletImageType > ChangeInformationFilterType;
		ChangeInformationFilterType::Pointer changeInformationFilter = ChangeInformationFilterType::New();


		std::cout<< "Projection of the higher level using Contour..\n";
		changeInformationFilter->SetInput( edgeContourFilter->GetOutput());

		const WaveletImageType::SpacingType spacing( 1 );
		changeInformationFilter->SetOutputSpacing( spacing);
		changeInformationFilter->ChangeSpacingOn();

		segmentationOnSource->SetNthInput(1, changeInformationFilter->GetOutput());

		
		RescaleIntensityImageFilterType::Pointer intensityImageFilter = RescaleIntensityImageFilterType::New();

		intensityImageFilter->SetInput(inputGradient->GetOutput());

		try
		{
			intensityImageFilter->Update();
		}
		catch (itk::ExceptionObject& e)
		{
			std::cerr << e << std::endl;
		}

		int maximumIntensity = intensityImageFilter->GetInputMaximum();

		char projectionEquation[100];
		
		segmentationOnSource->SetNthInput(1, changeInformationFilter->GetOutput());

		sprintf(projectionEquation, "if(b2 > 0, %f + b1, b1)", static_cast<float>(maximumIntensity * multiscaleContourRatioOfMaximumIntensity));

		segmentationOnSource->SetExpression(projectionEquation);

		sprintf(filenameBuffer,"Level%dWSegmentationOnImage.tif", decimationLevel);
	 	outputWriter->SetFileName(filenameBuffer);
		outputWriter->SetInput(segmentationOnSource->GetOutput());
		try
		{
			outputWriter->Update();
		}
		catch (itk::ExceptionObject& e)
		{
			std::cerr << e << std::endl;
		}

	}


	// ********************************************************************************************
	std::cout << "Computing watershed ..\n";

	typedef otb::WatershedSegmentationFilter<WaveletImageType, WaveletImageType> WatershedFilterType;

	WatershedFilterType::Pointer watershed = WatershedFilterType::New();
	watershed->SetLevel(0.3);
	watershed->SetThreshold(0.001);

	if(decimationLevel < 2)
	{
		watershed->SetInput(segmentationOnSource->GetOutput());
	}	
	else
	{
		watershed->SetInput(inputGradient->GetOutput());
	}

	// ********************************************************************************************
	// Produit une image de la sortie du traitement Watershed

	typedef itk::Functor::ScalarToRGBPixelFunctor<unsigned long> ColorMapFunctorType;
	typedef itk::UnaryFunctorImageFilter<WaveletImageType, RGBImageType, ColorMapFunctorType> ColorMapFilterType;
	ColorMapFilterType::Pointer colormapper = ColorMapFilterType::New();

	colormapper->SetInput(watershed->GetOutput());
							
	sprintf(filenameBuffer,"Level%dWatershed.tif", decimationLevel);
	postColorWatershedWriter->SetFileName(filenameBuffer);
	postColorWatershedWriter->SetInput(colormapper->GetOutput());
	try
	{
		postColorWatershedWriter->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}


	ImageToLabelMapFilterType::Pointer image2LabelMap =  ImageToLabelMapFilterType::New();
	
	image2LabelMap->SetInput(watershed->GetOutput());

	try
	{
		image2LabelMap->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}
	
  

	// ********************************************************************************************
	// Fusionne les petits objets

	typedef otb::ShapeAttributesLabelMapFilter<LabelMapType> ShapeLabelMapFilterType;

	ShapeLabelMapFilterType::Pointer shapeLabelMapFilter = ShapeLabelMapFilterType::New();
	shapeLabelMapFilter->SetInput(image2LabelMap->GetOutput());
	shapeLabelMapFilter->GetOutput()->SetAdjacencyMap(image2LabelMap->GetOutput()->GetAdjacencyMap());

	try
	{
		shapeLabelMapFilter->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}

	std::cout << "Object count before shape size filter : " << shapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << "\n";

	ShapeLabelMapFilterType::Pointer mergingShapeLabelMapFilter = shapeLabelMapFilter;

	for(int iterationCount = 0; iterationCount < mergeIterationCount; ++iterationCount)
	{
		std::cout << "Shape merge iteration " << iterationCount << "\n";

		auto labelObjectAdjacencyMap = mergingShapeLabelMapFilter->GetOutput()->GetAdjacencyMap();

		for(auto labelObjectIterator = labelObjectAdjacencyMap.begin(); 
				labelObjectIterator != labelObjectAdjacencyMap.end();  
				++labelObjectIterator)
		{
			int labelObjectMasterIndex = labelObjectIterator->first;
			float masterSizeValue = 0;

			try
			{
				masterSizeValue =  mergingShapeLabelMapFilter->GetOutput()->GetLabelObject(labelObjectMasterIndex)->GetAttribute("SHAPE::PhysicalSize");
		
			}
			catch (itk::ExceptionObject& e)
			{
				// Un "Throw" est toléré comme nous n'avons pas de contrôle sur l'implémentation de l'itérateur
				// qui apparaît comme incomplète.
				(void)e;
				// std::cerr << e << " for LabelObject " << labelObjectMasterIndex << std::endl;
			}

#if 1
// Cette algorithme n'est pas idéal mais produit des résultats satisfaisant dans le contexte,
// vue les problèmes d'itérateurs observés.

			for(auto LableObjectNeighborIterator = labelObjectIterator->second.begin(); 
				LableObjectNeighborIterator != labelObjectIterator->second.end(); 
				++LableObjectNeighborIterator)
			{
				int labelObjectNeighborIndex = *LableObjectNeighborIterator;

				float neighborSizeValue = 0;

				try
				{
					neighborSizeValue =  mergingShapeLabelMapFilter->GetOutput()->GetLabelObject(labelObjectNeighborIndex)->GetAttribute("SHAPE::PhysicalSize");
				}
				catch (itk::ExceptionObject& e)
				{
					// Un "Throw" est toléré comme nous n'avons pas de contrôle sur l'implémentation de l'itérateur
					// qui apparaît comme incomplète.
					(void)e;
					// std::cerr << e << std::endl;
				}

				 // Fusion des petits objets
				if(
					( neighborSizeValue < (labelObjectMinimumPixelSize * (std::pow(2.0, 2 - decimationLevel))))
					)
				{
					++LableObjectNeighborIterator;
					if(
						(masterSizeValue > (neighborSizeValue * 10 ) ) || // Try to merge to a big neighbor
						(LableObjectNeighborIterator == labelObjectIterator->second.end()) // Merge if this is the last neighbor
						) 
					{
						try
						{
							mergingShapeLabelMapFilter->GetOutput()->MergeLabels(labelObjectMasterIndex, labelObjectNeighborIndex);

						}
						catch (itk::ExceptionObject& e)
						{
							// Un "Throw" est toléré comme nous n'avons pas de contrôle sur l'implémentation de l'itérateur
							// qui apparaît comme incomplète.
							(void)e;
							// std::cerr << e << std::endl;
						}
					}
					--LableObjectNeighborIterator;
				}
			}
#else
// Pour une raison inconnu, cet algorithm ne fonctionne pas bien.
// C'est probablement lié au fait que l'implémentation des itérateurs semble incomplète.

			if(masterSizeValue < (labelObjectMinimumPixelSize * (std::pow(2.0, 2 - decimationLevel))))
			{
				// Find the biggest neighbor
				float biggestNeighborSize = 0.0f;
				int biggestNeighborIndex = -1;

				for(auto LableObjectNeighborIterator = labelObjectIterator->second.begin(); 
					LableObjectNeighborIterator != labelObjectIterator->second.end(); 
					++LableObjectNeighborIterator)
				{
					int labelObjectNeighborIndex = *LableObjectNeighborIterator;

					float neighborSizeValue = 0.0f;

					try
					{
						neighborSizeValue =  mergingShapeLabelMapFilter->GetOutput()->GetLabelObject(labelObjectNeighborIndex)->GetAttribute("SHAPE::PhysicalSize");

						if(neighborSizeValue > biggestNeighborSize)
						{
							biggestNeighborSize = neighborSizeValue;
							biggestNeighborIndex = labelObjectNeighborIndex;
						}
					}
					catch (itk::ExceptionObject& e)
					{
						(void)e;
						// std::cerr << e << std::endl;
					}
				}

				if(biggestNeighborIndex >= 0 )
				{
					try
					{
						mergingShapeLabelMapFilter->GetOutput()->MergeLabels(labelObjectMasterIndex, biggestNeighborIndex);
					}
					catch (itk::ExceptionObject& e)
					{
						(void)e;
						// std::cerr << e << std::endl;
					}
				}
			}
#endif
		}


		std::cout << "    Object count after shape size filter : " << mergingShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << "\n";
	}


	LabelMapToLabelImageFilterType::Pointer shapeLabelImage = LabelMapToLabelImageFilterType::New();
	shapeLabelImage->SetInput(mergingShapeLabelMapFilter->GetOutput());

	colormapper->SetInput(shapeLabelImage->GetOutput());

	sprintf(filenameBuffer,"Level%dObjectLabelShapeFiltered.tif", decimationLevel);
	postColorWatershedWriter->SetFileName(filenameBuffer);
	postColorWatershedWriter->SetInput(colormapper->GetOutput());
	try
	{
		postColorWatershedWriter->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}

	// ********************************************************************************************
	std::cout << "Computing texture indices ..\n";

	typedef otb::ScalarImageToTexturesFilter<WaveletImageType, WaveletImageType> ScalarImageToTexturesFilterType;

	ScalarImageToTexturesFilterType::Pointer scalarImageToTexturesFilter = ScalarImageToTexturesFilterType::New();

	scalarImageToTexturesFilter->SetInput(inputReader->GetOutput());

	typedef WaveletImageType::SizeType RadiusSizeType;
	RadiusSizeType sradius;
	sradius.Fill(textureRadius);

	scalarImageToTexturesFilter->SetRadius(sradius);

	typedef WaveletImageType::OffsetType OffsetType;
	OffsetType offset;
	offset[0] =  0;
	offset[1] =  0;
	scalarImageToTexturesFilter->SetOffset(offset);


	RescaleIntensityImageFilterType::Pointer intensityTextureImageFilter = RescaleIntensityImageFilterType::New();

	intensityTextureImageFilter->SetInput(inputReader->GetOutput());

	try
	{
		intensityTextureImageFilter->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}

	float maximumIntensity = intensityTextureImageFilter->GetInputMaximum();

	scalarImageToTexturesFilter->SetInputImageMinimum(0);
    scalarImageToTexturesFilter->SetInputImageMaximum(maximumIntensity);

	try
	{
		scalarImageToTexturesFilter->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}

	// *** Produit l'image de l'indice "Cluster Shade" ***

	intensityTextureImageFilter->SetInput(scalarImageToTexturesFilter->GetClusterShadeOutput());

	sprintf(filenameBuffer,"Level%dIndiceClusterShade.tif", decimationLevel);
	outputWriter->SetFileName(filenameBuffer);
	outputWriter->SetInput(intensityTextureImageFilter->GetOutput());
	try
	{
		outputWriter->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}
	
	float indiceClusterShadeMaximum = intensityTextureImageFilter->GetInputMaximum();

	// *** Produit l'image de l'indice "Energy" ***

	intensityTextureImageFilter->SetInput(scalarImageToTexturesFilter->GetEnergyOutput());

	sprintf(filenameBuffer,"Level%dIndiceEnergy.tif", decimationLevel);
	outputWriter->SetFileName(filenameBuffer);
	outputWriter->SetInput(intensityTextureImageFilter->GetOutput());
	try
	{
		outputWriter->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}

	float indiceEnergyMaximum = intensityTextureImageFilter->GetInputMaximum();

	// ********************************************************************************************
	std::cout << "Computing advanced texture indices ..\n";

	typedef otb::ScalarImageToAdvancedTexturesFilter<WaveletImageType, WaveletImageType> ScalarImageToAdvancedTexturesFilterType;

	ScalarImageToAdvancedTexturesFilterType::Pointer scalarImageToAdvancedTexturesFilter = ScalarImageToAdvancedTexturesFilterType::New();

	scalarImageToAdvancedTexturesFilter->SetInput(inputReader->GetOutput());

	sradius.Fill(textureRadius);

	scalarImageToAdvancedTexturesFilter->SetRadius(sradius);

	scalarImageToAdvancedTexturesFilter->SetOffset(offset);

	intensityTextureImageFilter->SetInput(inputReader->GetOutput());

	try
	{
		intensityTextureImageFilter->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}

	maximumIntensity = intensityTextureImageFilter->GetInputMaximum();

	scalarImageToAdvancedTexturesFilter->SetInputImageMinimum(0);
    scalarImageToAdvancedTexturesFilter->SetInputImageMaximum(maximumIntensity);

	try
	{
		scalarImageToAdvancedTexturesFilter->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}


	// *** Produit l'image de l'indice "Mean" ***

	intensityTextureImageFilter->SetInput(scalarImageToAdvancedTexturesFilter->GetMeanOutput());

	sprintf(filenameBuffer,"Level%dIndiceMean.tif", decimationLevel);
	outputWriter->SetFileName(filenameBuffer);
	outputWriter->SetInput(intensityTextureImageFilter->GetOutput());
	try
	{
		outputWriter->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}

	float indiceMeanMaximum = intensityTextureImageFilter->GetInputMaximum();
	

	// *** Produit l'image de l'indice "Sum Variance" ***

	intensityTextureImageFilter->SetInput(scalarImageToAdvancedTexturesFilter->GetSumVarianceOutput());

	sprintf(filenameBuffer,"Level%dIndiceSumVariance.tif", decimationLevel);
	outputWriter->SetFileName(filenameBuffer);
	outputWriter->SetInput(intensityTextureImageFilter->GetOutput());
	try
	{
		outputWriter->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}

	float indiceSumVarianceMaximum = intensityTextureImageFilter->GetInputMaximum();


	// ********************************************************************************************
	
	StatisticsLabelMapFilterType::Pointer statisticLabelMapFilter = StatisticsLabelMapFilterType::New();

	// ****************************************************************************************************************


	for(unsigned int currentMergeIteration = 0; currentMergeIteration < mergeIterationCount; currentMergeIteration++)
	{
		std::cout << "Merging Iteration " << currentMergeIteration << "\n";

		statisticLabelMapFilter->SetInput(mergingShapeLabelMapFilter->GetOutput());

		// ********************************************************************************************
		// Calcul les moyennes d'indice de ClusterShade.

		statisticLabelMapFilter->SetFeatureImage(scalarImageToTexturesFilter->GetClusterShadeOutput());

		try
		{
			statisticLabelMapFilter->Update();
		}
		catch (itk::ExceptionObject& e)
		{
			std::cerr << e << std::endl;
		}


		auto labelObjectAdjacencyMap = mergingShapeLabelMapFilter->GetOutput()->GetAdjacencyMap();

		// Renomme l'attribut
		for(auto labelObjectIterator = labelObjectAdjacencyMap.begin(); 
			labelObjectIterator != labelObjectAdjacencyMap.end();  
			++labelObjectIterator)
		{
			int labelObjectMasterIndex = labelObjectIterator->first;	

			try
			{
				float masterValue = statisticLabelMapFilter->GetOutput()->GetLabelObject(labelObjectMasterIndex)->GetAttribute("STATS::Default::Mean");

				mergingShapeLabelMapFilter->GetOutput()->GetLabelObject(labelObjectMasterIndex)->SetAttribute("Multiechelle::ClusterShade", masterValue / indiceClusterShadeMaximum);
			}
			catch (itk::ExceptionObject& e)
			{
				std::cerr << e << std::endl;
			}
		}

		// ********************************************************************************************
		// Calcul les moyennes d'indice de Energy.

		statisticLabelMapFilter->SetFeatureImage(scalarImageToTexturesFilter->GetEnergyOutput());

		try
		{
			statisticLabelMapFilter->Update();
		}
		catch (itk::ExceptionObject& e)
		{
			std::cerr << e << std::endl;
		}


		labelObjectAdjacencyMap = mergingShapeLabelMapFilter->GetOutput()->GetAdjacencyMap();

		// Renomme l'attribut
		for(auto labelObjectIterator = labelObjectAdjacencyMap.begin(); 
			labelObjectIterator != labelObjectAdjacencyMap.end();  
			++labelObjectIterator)
		{
			int labelObjectMasterIndex = labelObjectIterator->first;	

			try
			{
				float masterValue = statisticLabelMapFilter->GetOutput()->GetLabelObject(labelObjectMasterIndex)->GetAttribute("STATS::Default::Mean");

				mergingShapeLabelMapFilter->GetOutput()->GetLabelObject(labelObjectMasterIndex)->SetAttribute("Multiechelle::Energy", masterValue / indiceEnergyMaximum);
			}
			catch (itk::ExceptionObject& e)
			{
				std::cerr << e << std::endl;
			}
		}
		
		// ********************************************************************************************
		// Calcul les moyennes d'indice de Mean.

		statisticLabelMapFilter->SetFeatureImage(scalarImageToAdvancedTexturesFilter->GetMeanOutput());

		try
		{
			statisticLabelMapFilter->Update();
		}
		catch (itk::ExceptionObject& e)
		{
			std::cerr << e << std::endl;
		}


		labelObjectAdjacencyMap = mergingShapeLabelMapFilter->GetOutput()->GetAdjacencyMap();

		// Renomme l'attribute
		for(auto labelObjectIterator = labelObjectAdjacencyMap.begin(); 
			labelObjectIterator != labelObjectAdjacencyMap.end();  
			++labelObjectIterator)
		{
			int labelObjectMasterIndex = labelObjectIterator->first;	

			try
			{
				float masterValue = statisticLabelMapFilter->GetOutput()->GetLabelObject(labelObjectMasterIndex)->GetAttribute("STATS::Default::Mean");

				mergingShapeLabelMapFilter->GetOutput()->GetLabelObject(labelObjectMasterIndex)->SetAttribute("Multiechelle::Mean", masterValue / indiceMeanMaximum);
			}
			catch (itk::ExceptionObject& e)
			{
				std::cerr << e << std::endl;
			}
		}

		// ********************************************************************************************
		// Calcul les moyennes d'indice de SumVariance.

		statisticLabelMapFilter->SetFeatureImage(scalarImageToAdvancedTexturesFilter->GetSumVarianceOutput());

		try
		{
			statisticLabelMapFilter->Update();
		}
		catch (itk::ExceptionObject& e)
		{
			std::cerr << e << std::endl;
		}


		labelObjectAdjacencyMap = mergingShapeLabelMapFilter->GetOutput()->GetAdjacencyMap();

		// Renomme l'attribut.
		for(auto labelObjectIterator = labelObjectAdjacencyMap.begin(); 
			labelObjectIterator != labelObjectAdjacencyMap.end();  
			++labelObjectIterator)
		{
			int labelObjectMasterIndex = labelObjectIterator->first;	

			try
			{
				float masterValue = statisticLabelMapFilter->GetOutput()->GetLabelObject(labelObjectMasterIndex)->GetAttribute("STATS::Default::Mean");

				mergingShapeLabelMapFilter->GetOutput()->GetLabelObject(labelObjectMasterIndex)->SetAttribute("Multiechelle::SumVariance", masterValue / indiceSumVarianceMaximum);
			}
			catch (itk::ExceptionObject& e)
			{
				std::cerr << e << std::endl;
			}
		}
		// ********************************************************************************************
		// Effecture la fusion basé sur les indices de textures.

		// std::cout << "    LabelObject count before merging : " << mergingShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << "\n";
		// PrintLabelObjectAttributes(mergingShapeLabelMapFilter->GetOutput()->GetNthLabelObject(1));

		labelObjectAdjacencyMap = mergingShapeLabelMapFilter->GetOutput()->GetAdjacencyMap();

		for(auto labelObjectIterator = labelObjectAdjacencyMap.begin(); 
			labelObjectIterator != labelObjectAdjacencyMap.end();  
			++labelObjectIterator)
		{

			int labelObjectMasterIndex = labelObjectIterator->first;
			float masterEnergyValue = 0;
			float masterClusterShadeValue = 0;
			float masterMeanValue = 0;
			float masterSumVarianceValue = 0;

			try
			{
				masterEnergyValue = mergingShapeLabelMapFilter->GetOutput()->GetLabelObject(labelObjectMasterIndex)->GetAttribute("Multiechelle::Energy");
				masterClusterShadeValue = mergingShapeLabelMapFilter->GetOutput()->GetLabelObject(labelObjectMasterIndex)->GetAttribute("Multiechelle::ClusterShade");
				masterMeanValue = mergingShapeLabelMapFilter->GetOutput()->GetLabelObject(labelObjectMasterIndex)->GetAttribute("Multiechelle::Mean");
				masterSumVarianceValue = mergingShapeLabelMapFilter->GetOutput()->GetLabelObject(labelObjectMasterIndex)->GetAttribute("Multiechelle::SumVariance");

			}
			catch (itk::ExceptionObject& e)
			{
				// Un "Throw" est toléré comme nous n'avons pas de contrôle sur l'implémentation de l'itérateur
				// qui apparaît comme incomplète.
				(void)e;
				// std::cerr << e << " for LabelObject " << labelObjectMasterIndex << std::endl;
			}

			for(auto LableObjectNeighborIterator = labelObjectIterator->second.begin(); 
				LableObjectNeighborIterator != labelObjectIterator->second.end(); 
				++LableObjectNeighborIterator)
			{
				int labelObjectNeighborIndex = *LableObjectNeighborIterator;

				float neighborEnergyValue = 0;
				float neighborClusterShadeValue = 0;
				float neighborMeanValue = 0;
				float neighborSumVarianceValue = 0;

				try
				{
					neighborEnergyValue = mergingShapeLabelMapFilter->GetOutput()->GetLabelObject(labelObjectNeighborIndex)->GetAttribute("Multiechelle::Energy");
					neighborClusterShadeValue = mergingShapeLabelMapFilter->GetOutput()->GetLabelObject(labelObjectNeighborIndex)->GetAttribute("Multiechelle::ClusterShade");
					neighborMeanValue = mergingShapeLabelMapFilter->GetOutput()->GetLabelObject(labelObjectMasterIndex)->GetAttribute("Multiechelle::Mean");
					neighborSumVarianceValue = mergingShapeLabelMapFilter->GetOutput()->GetLabelObject(labelObjectMasterIndex)->GetAttribute("Multiechelle::SumVariance");
				}
				catch (itk::ExceptionObject& e)
				{
					// Un "Throw" est toléré comme nous n'avons pas de contrôle sur l'implémentation de l'itérateur
					// qui apparaît comme incomplète.
					(void)e;
					// std::cerr << e << std::endl;
				}


				float masterValue =  masterEnergyValue + masterClusterShadeValue + masterMeanValue + masterSumVarianceValue;
					
				float neighborValue = neighborEnergyValue + neighborClusterShadeValue + neighborMeanValue + neighborSumVarianceValue;

				float currentIndiceDelta = std::fabs(masterValue - neighborValue);

				if( currentIndiceDelta < mergingThreshold) 
				{
					try
					{
						mergingShapeLabelMapFilter->GetOutput()->MergeLabels(labelObjectMasterIndex, labelObjectNeighborIndex);

					}
					catch (itk::ExceptionObject& e)
					{
						// Un "Throw" est toléré comme nous n'avons pas de contrôle sur l'implémentation de l'itérateur
						// qui apparaît comme incomplète.
						(void)e;
						// std::cerr << e << std::endl;
					}
				}
			}
		}
		std::cout << "    Object count after merging : " << mergingShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << "\n";
	}
	
	

	// ****************************************************************************************************************
	// Produit la carte de segmentation.
	std::cout << "Extracting segmentation..\n";

	LabelMapToLabelImageFilterType::Pointer finalLabelImage = LabelMapToLabelImageFilterType::New();
	finalLabelImage->SetInput(mergingShapeLabelMapFilter->GetOutput());

	outputWriter->SetInput(finalLabelImage->GetOutput());////////////////////////////////////////////////////////////////////////////////////////////////////////////////////je pense que c'est ici qu'il faut agir
	outputWriter->SetFileName(outputFileName);
	try
	{
		outputWriter->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}

	colormapper->SetInput(finalLabelImage->GetOutput());

	sprintf(filenameBuffer,"Level%dObjectLabelMerged.tif", decimationLevel);
	postColorWatershedWriter->SetFileName(filenameBuffer);
	postColorWatershedWriter->SetInput(colormapper->GetOutput());
	try
	{
		postColorWatershedWriter->Update();
	}
	catch (itk::ExceptionObject& e)
	{
		std::cerr << e << std::endl;
	}

	// ****************************************************************************************************************
	// Supperpose la segmentation produite sur l'image d'origine si requis.

	if(decimationLevel == 0)
	{
		std::cout << "Creating segmented image..\n";

		LabelImageContourFilterType::Pointer edgeFilter = LabelImageContourFilterType::New();
		edgeFilter->SetInput(finalLabelImage->GetOutput());
		edgeFilter->SetFullyConnected(false);

		BinaryThinningImageFilterType::Pointer edgeThinnerFilter = BinaryThinningImageFilterType::New();
		edgeThinnerFilter->SetInput(edgeFilter->GetOutput());

		typedef otb::BandMathImageFilter<WaveletImageType> BandMathFilterType;
		BandMathFilterType::Pointer segmentationOnSource = BandMathFilterType::New();

		segmentationOnSource->SetNthInput(0, inputReader->GetOutput());
		segmentationOnSource->SetNthInput(1, edgeThinnerFilter->GetOutput());


		RescaleIntensityImageFilterType::Pointer intensityImageFilter = RescaleIntensityImageFilterType::New();

		intensityImageFilter->SetInput(inputReader->GetOutput());

		try
		{
			intensityImageFilter->Update();
		}
		catch (itk::ExceptionObject& e)
		{
			std::cerr << e << std::endl;
		}

		int maximumIntensity = intensityImageFilter->GetInputMaximum();

		char projectionEquation[100];

		sprintf(projectionEquation, "if(b2 > 0, %f, b1)", static_cast<float>(maximumIntensity*1.3f));

		segmentationOnSource->SetExpression(projectionEquation);

	 	outputWriter->SetFileName("Segmentation.tif");
		outputWriter->SetInput(segmentationOnSource->GetOutput());
		try
		{
			outputWriter->Update();
		}
		catch (itk::ExceptionObject& e)
		{
			std::cerr << e << std::endl;
		}


		segmentationOnSource->Update();
		RescaleIntensityImageFilterType1::Pointer intensityImageFilter1 = RescaleIntensityImageFilterType1::New();
        intensityImageFilter1->SetInput(segmentationOnSource->GetOutput());
		intensityImageFilter1->Update();

		SetParameterOutputImage("out",intensityImageFilter1->GetOutput());
		//otb::Wrapper::Application::SetParameterOutputImagePixelType
	}
	

	std::cout << "Segmentation DONE\n\n\n";
	return EXIT_SUCCESS;
}


     };

}
   
}
OTB_APPLICATION_EXPORT(otb::Wrapper::MultiscaleSegmentation);