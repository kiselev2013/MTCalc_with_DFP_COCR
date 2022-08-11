/**
 *  This code is freely available under the following conditions:
 *
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *
 *
 * The file contains
 * 1. Regular parallelepiped finite element mesh building functions using given settings
 * 2. Functions for writing mesh files, settings and task parameters files needed for calculation
 *
 *  Written by Ph.D. Dmitry S. Kiselev
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  harlequin_00@mail.ru
 *  Version 1.1 July, 2021
*/

#include <stdio.h>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include "BasicOperations.h"
#include <chrono>

// Program start time
std::chrono::high_resolution_clock::time_point tProgramStart;

// Gate for determining if materials properties are equal or not
#define MaterialEqualEps 1e-7

// Types of material properties values
enum PhysicalValueType
{
	Constant,
	Tensor,
	Vector,
	Table,
	Curve
};

// Names of material properties
enum PhysicalValueName
{
	SigmaN,
	SigmaH,
	SigmaV,
	Eps,
	Mu
};

// Material property class
class PhysicalValue
{
public:
	PhysicalValueType valueType;
	double constantValue;

	PhysicalValue::PhysicalValue()
	{
		valueType = PhysicalValueType::Constant;
		constantValue = 0.0;
	}
	PhysicalValue::PhysicalValue(const PhysicalValue &value)
	{
		valueType = value.valueType;
		constantValue = value.constantValue;
	}

	// Check if two constant values are equal
	bool ConstantEqual(double c1, double c2)
	{
		if (fabs(c1) < 1e-30 && fabs(c2) < 1e-30)
			return true;
		if (fabs(1.0 - min(c2, c1) / max(c2, c1)) < MaterialEqualEps)
			return true;
		return false;
	}

	// Check if two property value and this are equal
	bool PhysicalValue::Equal(PhysicalValue &value, bool compareTensors = false)
	{
		if (value.valueType != valueType)
			return false;

		switch (valueType)
		{
		case PhysicalValueType::Constant:
			return ConstantEqual(constantValue, value.constantValue);
		case PhysicalValueType::Vector:
			break;
		case PhysicalValueType::Tensor:
			break;
		case PhysicalValueType::Table:
			break;
		case PhysicalValueType::Curve:
			break;
		}

		return false;
	}

};

// Material class
class Material
{
public:

	// Material properties:
	PhysicalValue SigmaH;
	PhysicalValue SigmaV;
	PhysicalValue SigmaN;
	PhysicalValue Eps;
	PhysicalValue Mu;

	Material::Material()
	{
		Mu.constantValue = 1.0;
		Eps.constantValue = 0.0;
	}
	Material::Material(const Material &material)
	{
		CopyPhysicals(material);
	}

	// Copy material properties from another material to this
	void Material::CopyPhysicals(const Material &material)
	{
		this->SigmaH = material.SigmaH;
		this->SigmaV = material.SigmaV;
		this->SigmaN = material.SigmaN;
		this->Eps = material.Eps;
		this->Mu = material.Mu;
	}

	// Check if another material is equal to this
	bool Material::Equal(Material &material)
	{
		if (!material.SigmaH.Equal(SigmaH)) return false;
		if (!material.SigmaV.Equal(SigmaV)) return false;
		if (!material.SigmaN.Equal(SigmaN)) return false;
		if (!material.Eps.Equal(Eps)) return false;
		if (!material.Mu.Equal(Mu)) return false;

		return true;
	}

	PhysicalValue& Material::GetValue(PhysicalValueName name)
	{
		switch (name)
		{
		case PhysicalValueName::SigmaN:
			return SigmaN;
		case PhysicalValueName::SigmaH:
			return SigmaH;
		case PhysicalValueName::SigmaV:
			return SigmaV;
		case PhysicalValueName::Eps:
			return Eps;
		case PhysicalValueName::Mu:
			return Mu;
		}
	}
};

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));

}
using namespace std;

// Common settings class
class Settings
{
public:
	double steps[3]; // Mesh initial steps
	double sparse[3]; // Sparce coefficient
	double farBound[3]; // Distance from receivers bounding box+gap to the boundary of the calculation domain
	double gapFromReceivers[3]; // Gap between receivers bounding box and sparced mesh area
	double eps = 0.1; // Gate for determining if coordinate lines are equal
	int currentDirection = 1; // Current direction for current calculation
	double frequency; // Frequency for current calculation
	int threadsCount = 1;  // Threads count for parallel calculation
	double SLAESolutionEps = 1e-4;  // relative residual for determining if SLAE solving is finished
	int maxIterationsCount = 10000; // Maximum iterations when solving SLAE

	Settings() { }
};

// Box class for different purposes
template <typename T>
class GeoBox
{
public:

	T coordinates[6]; // x0, x1, y0, y1, z0, z1
	int material = 0; // Material number
	Material materialValue; // Matrial value

	GeoBox()
	{

	}
};

// Layer class
template <typename T>
class GeoLayer
{
public:

	T top; // Upper layer coordinate
	T bottom; // Lower layer coordinate
	int material = 0;  // Material number
	Material materialValue; // Material value

	GeoLayer()
	{

	}
};

// Receiver class
class Receiver
{
public:

	double coordinates[3]; // x, y, z

	Receiver()
	{

	}
};


// Find value in v
int BinarySearchInVectorSorted(vector<double> &v, double value)
{
	int b[2] = { 0, v.size() - 1 };
	int m;

	if (value < v.front()) return 0;
	if (value > v.back()) return v.size()-1;

	while (b[1] - b[0] > 1)
	{
		m = (b[1] + b[0]) / 2;
		if (value > v[m])
			b[0] = m;
		else
			b[1] = m;
	}

	if (b[0] == b[1])
		return b[0];

	if (fabs(v[b[0]] - value) < fabs(v[b[1]] - value))
		return b[0];

	return b[1];
}

// Build 1D mesh within bounds with constant step
int Build1DUniformMesh(double step, double *bounds, vector<double> &mesh)
{
	try
	{
		if (bounds[1] < bounds[0]) { write_to_log("Error : Build1DUniformMesh : Mesh building bounds are inadequate\n"); return 1; }
		if (step < 0) { write_to_log("Error : Build1DUniformMesh : Step is less tha zero\n"); return 1; }
		int stepsCount = (bounds[1] - bounds[0]) / step + 0.5;

		if (fabs(step * stepsCount - (bounds[1] - bounds[0])) > 1e-5)
			stepsCount++;

		step = (bounds[1] - bounds[0]) / stepsCount;

		mesh.resize(stepsCount + 1);
		for (size_t i = 0; i < stepsCount + 1; i++)
			mesh[i] = bounds[0] + i*step;

		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : Build1DUniformMesh : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Build 1D mesh within bounds with increasing by sparce step. Direction > 0 is from left to right, otherwise - from right to left
int Build1DSparsedMesh(double step, double sparse, int direction, double *bounds, vector<double> &mesh)
{
	try
	{
		if (bounds[1] < bounds[0]) { write_to_log("Error : Build1DSparsedMesh : Mesh building bounds are inadequate\n"); return 1; }
		if (step < 0) { write_to_log("Error : Build1DSparsedMesh : Step is less tha zero\n"); return 1; }
		if (sparse < 1.0) { write_to_log("Error : Build1DSparsedMesh : Sparse coefficient is less tha 1.0\n"); return 1; }

		double coord;
		if (direction > 0)
		{
			coord = bounds[0];
			mesh.push_back(coord);
			do {
				coord += step;
				step *= sparse;
				mesh.push_back(coord);
			} while (coord < bounds[1]);
		}
		else
		{
			vector<double> meshReverse;
			coord = bounds[1];
			meshReverse.push_back(coord);
			do {
				coord -= step;
				step *= sparse;
				meshReverse.push_back(coord);
			} while (coord > bounds[0]);

			mesh.resize(meshReverse.size());
			int meshIndex = 0;
			for (int meshReverseIndex = meshReverse.size() - 1; meshReverseIndex >= 0; meshReverseIndex--, meshIndex++)
				mesh[meshIndex] = meshReverse[meshReverseIndex];
		}

		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : Build1DSparsedMesh : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Build base 1d mesh using model bounds
int BuildBaseMesh1D(double step, double sparse, double *modelBounds, vector<double> &mesh)
{
	try
	{
		vector<double> meshL, meshM, meshR;

		if (Build1DUniformMesh(step, modelBounds + 1, meshM) != 0) { write_to_log("Error : BuildBaseMesh1D : Could not build middle uniform mesh\n"); return 1; }
		if (Build1DSparsedMesh(step, sparse, -1, modelBounds + 0, meshL) != 0) { write_to_log("Error : BuildBaseMesh1D : Could not build left mesh\n"); return 1; }
		if (Build1DSparsedMesh(step, sparse,  1, modelBounds + 2, meshR) != 0) { write_to_log("Error : BuildBaseMesh1D : Could not build roght mesh\n"); return 1; }

		mesh.resize(meshL.size() + meshM.size() + meshR.size() - 2);

		int meshIndex = 0;
		for (int tmpIndex = 0; tmpIndex < meshL.size() - 1; tmpIndex++, meshIndex++)
			mesh[meshIndex] = meshL[tmpIndex];

		for (int tmpIndex = 0; tmpIndex < meshM.size() - 1; tmpIndex++, meshIndex++)
			mesh[meshIndex] = meshM[tmpIndex];

		for (int tmpIndex = 0; tmpIndex < meshR.size(); tmpIndex++, meshIndex++)
			mesh[meshIndex] = meshR[tmpIndex];

		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : BuildBaseMesh1D : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Build base 3d mesh using model bounds
int BuildBaseMesh3D(Settings &meshSettings, double **modelBounds, vector<double> *mesh)
{
	try
	{
		int direction;
		direction = 0; if (BuildBaseMesh1D(meshSettings.steps[direction], meshSettings.sparse[direction], modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : BuildBaseMesh3D : Could not build X base mesh\n"); return 1; }
		direction = 1; if (BuildBaseMesh1D(meshSettings.steps[direction], meshSettings.sparse[direction], modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : BuildBaseMesh3D : Could not build Y base mesh\n"); return 1; }
		direction = 2; if (BuildBaseMesh1D(meshSettings.steps[direction], meshSettings.sparse[direction], modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : BuildBaseMesh3D : Could not build Z base mesh\n"); return 1; }
		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : BuildBaseMesh3D : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Build template base 3d mesh using model bounds
int BuildBaseMesh3DTemplate(double **modelBounds, vector<double> *mesh, vector<double> *meshTemplate)
{
	try
	{
		
		for (int direction = 0; direction < 3; direction++)
		{
			
			auto &mesh1D = mesh[direction];
			auto &mesh1DTemplate= meshTemplate[direction];
			mesh1DTemplate.resize(mesh1D.size());
			modelBounds[direction][0] = mesh1D.front();
			modelBounds[direction][3] = mesh1D.back();

			for (int meshIndex = 0; meshIndex < mesh1D.size(); meshIndex++)
			{
				double coordinate = mesh1D[meshIndex];
				if (coordinate < modelBounds[direction][1])
				{
					mesh1DTemplate[meshIndex] = max(0.0, (coordinate - modelBounds[direction][0]) / (modelBounds[direction][1] - modelBounds[direction][0]));
				}
				else if (coordinate < modelBounds[direction][2])
				{
					mesh1DTemplate[meshIndex] = (coordinate - modelBounds[direction][1]) / (modelBounds[direction][2] - modelBounds[direction][1]) + 1;
				}
				else
					mesh1DTemplate[meshIndex] = min(1.0, (coordinate - modelBounds[direction][2]) / (modelBounds[direction][3] - modelBounds[direction][2])) + 2;
			}
		}

		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : BuildBaseMesh3DTemplate : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Align 1D mesh begin and end with bounds 
int AlignMeshBounds(double *bounds, vector<double> &mesh)
{
	try
	{
		mesh[0] = bounds[0];
		if ((mesh[1] - mesh[0]) < (mesh[2] - mesh[1]) * 0.3)
			mesh.erase(mesh.begin() + 1);

		mesh[mesh.size() - 1] = bounds[3];
		if ((mesh[mesh.size() - 1] - mesh[mesh.size() - 2]) < (mesh[mesh.size() - 2] - mesh[mesh.size() - 3]) * 0.3)
			mesh.erase(mesh.begin() + mesh.size() - 2);

		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : AlignMeshBounds : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Align 3D mesh begin and end with bounds 
int AlignMeshBounds(double **modelBounds, vector<double> *mesh)
{
	try
	{
		int direction;
		direction = 0; if (AlignMeshBounds(modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : AlignMeshBounds : Could not align X mesh\n"); return 1; }
		direction = 1; if (AlignMeshBounds(modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : AlignMeshBounds : Could not align Y mesh\n"); return 1; }
		direction = 2; if (AlignMeshBounds(modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : AlignMeshBounds : Could not align Z mesh\n"); return 1; }
		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : AlignMeshBounds : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Sort coordinates
int InsertCoordinateIntoSortedArray(vector<double> &arr, double value, double eps)
{
	try
	{
		if (value < arr.front())
			return 0;
		if (value > arr.back())
			return 0;
		int	foundIndex = BinarySearchInVectorSorted(arr, value);
		int neighbourIndex = arr[foundIndex] > value && foundIndex != 0 ? foundIndex - 1 : foundIndex + 1;
		double step = fabs(arr[foundIndex] - arr[neighbourIndex]);
		if (fabs(arr[foundIndex] - value) / step > eps)
		{
			arr.insert(arr.begin() + (arr[foundIndex] > value ? foundIndex : foundIndex + 1), value);
		}
		else
			arr[foundIndex] = value;

		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : InsertCoordinateIntoSortedArray : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Create additional mesh lines for objects
int InsertObjectsIntoMesh(vector<GeoBox<double>> &objects, double eps, vector<double> *mesh)
{
	try
	{
		char buf[2048];
		int objectNumber = 1;
		for (auto obj = objects.begin(); obj != objects.end(); ++obj, objectNumber++)
		{
			for (int bound = 0; bound < 6; bound++)
			{
				if (InsertCoordinateIntoSortedArray(mesh[bound / 2], obj->coordinates[bound], eps) != 0) 
				{
					sprintf(buf, "Error : InsertObjectsIntoMesh : Could not insert object %d bound %d into mesh\n", objectNumber, bound+1); write_to_log(buf);
					return 1; 
				}
			}
		}
		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : InsertObjectsIntoMesh : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Create additional mesh lines for layers
int InsertLayersIntoMesh(vector<GeoLayer<double>> &layers, double eps, vector<double> *mesh)
{
	try
	{
		char buf[2048];
		int layerNumber = 1;
		for (auto layer = layers.begin(); layer != layers.end(); ++layer, layerNumber++)
		{
			if (layerNumber == layers.size())
				break;
			if (InsertCoordinateIntoSortedArray(mesh[2], layer->bottom, eps) != 0)
			{
				sprintf(buf, "Error : InsertLayersIntoMesh : Could not insert object %d bound %d into mesh\n", layerNumber, layer->bottom); write_to_log(buf);
				return 1;
			}
		}
		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : InsertLayersIntoMesh : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Round double value with decrement
double RoundDouble_(double value, double decrement)
{
	int n = fabs(value / decrement);
	return value < 0 ? -(n + 1) * decrement : n * decrement;
}

// Calculate receivers bounding box
int CalculateReceiversBounds(vector<Receiver> &receivers, Settings &meshSettings, GeoBox<double> &receiversBounds)
{
	try
	{
		if (receivers.size() == 0) { write_to_log("Error : CalculateReceiversBounds : No receivers\n"); return 1; }

		receiversBounds.coordinates[0] = receiversBounds.coordinates[2] = receiversBounds.coordinates[4] = DBL_MAX;
		receiversBounds.coordinates[1] = receiversBounds.coordinates[3] = receiversBounds.coordinates[5] = -DBL_MAX;

		for (auto receiver = receivers.begin(); receiver != receivers.end(); ++receiver)
		{
			for (int direction = 0; direction < 3; direction++)
			{
				if (receiversBounds.coordinates[direction * 2 + 0] > receiver->coordinates[direction]) receiversBounds.coordinates[direction * 2 + 0] = receiver->coordinates[direction];
				if (receiversBounds.coordinates[direction * 2 + 1] < receiver->coordinates[direction]) receiversBounds.coordinates[direction * 2 + 1] = receiver->coordinates[direction];
			}
		}

		for (int direction = 0; direction < 3; direction++)
		{
			receiversBounds.coordinates[direction * 2 + 0] -= meshSettings.gapFromReceivers[direction];
			receiversBounds.coordinates[direction * 2 + 1] += meshSettings.gapFromReceivers[direction];

			double d = receiversBounds.coordinates[direction * 2 + 1] - receiversBounds.coordinates[direction * 2 + 0];
			int stepsCount = (int)(d / meshSettings.steps[direction] + 0.5);
			if (fabs(stepsCount * meshSettings.steps[direction] - d) > 1e-7)
			{
				stepsCount++;
				receiversBounds.coordinates[direction * 2 + 0] = RoundDouble_(receiversBounds.coordinates[direction * 2 + 0], meshSettings.steps[direction]);
				receiversBounds.coordinates[direction * 2 + 1] = receiversBounds.coordinates[direction * 2 + 0] + stepsCount * meshSettings.steps[direction];
			}
		}

		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : CalculateReceiversBounds : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}

}

// Calculate receivers bounding box with rounding
int CalculateReceiversBounds(vector<Receiver> &receivers, GeoBox<double> &receiversBounds)
{
	try
	{
		if (receivers.size() == 0) { write_to_log("Error : CalculateReceiversBounds : No receivers\n"); return 1; }

		receiversBounds.coordinates[0] = receiversBounds.coordinates[2] = receiversBounds.coordinates[4] = DBL_MAX;
		receiversBounds.coordinates[1] = receiversBounds.coordinates[3] = receiversBounds.coordinates[5] = -DBL_MAX;

		for (auto receiver = receivers.begin(); receiver != receivers.end(); ++receiver)
		{
			for (int direction = 0; direction < 3; direction++)
			{
				if (receiversBounds.coordinates[direction * 2 + 0] > receiver->coordinates[direction]) receiversBounds.coordinates[direction * 2 + 0] = receiver->coordinates[direction];
				if (receiversBounds.coordinates[direction * 2 + 1] < receiver->coordinates[direction]) receiversBounds.coordinates[direction * 2 + 1] = receiver->coordinates[direction];
			}
		}

		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : CalculateReceiversBounds : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}

}

// Calculatemodel bounding box
int CalculateModelBounds(vector<Receiver> &receivers, Settings &meshSettings, double **bounds)
{
	try
	{
		GeoBox<double> receiversBoundsNotRounded;
		GeoBox<double> receiversBounds;
		if (CalculateReceiversBounds(receivers, receiversBoundsNotRounded) != 0) { write_to_log("Error : CalculateModelBounds : Could not calculate receiver bounds\n"); return 1; }
		if (CalculateReceiversBounds(receivers, meshSettings, receiversBounds) != 0) { write_to_log("Error : CalculateModelBounds : Could not calculate receiver bounds\n"); return 1; }

		for (int direction = 0; direction < 3; direction++)
		{
			bounds[direction][0] = receiversBoundsNotRounded.coordinates[direction * 2] - meshSettings.farBound[direction];
			bounds[direction][1] = receiversBounds.coordinates[direction * 2];
			bounds[direction][2] = receiversBounds.coordinates[direction * 2 + 1];
			bounds[direction][3] = receiversBoundsNotRounded.coordinates[direction * 2 + 1] + meshSettings.farBound[direction];
		}
		

		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : CalculateModelBounds : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Apply layers materials to mesh
int SetLayersMaterials(vector<GeoLayer<double>> &layers, vector<Material> &layerMaterials, vector<double> *mesh)
{
	try
	{
		int topIndex, bottomIndex, cell;
		int nx = mesh[0].size() - 1;
		int nxy = nx * (mesh[1].size() - 1);
		layerMaterials.resize(nxy * (mesh[2].size() - 1));

		for (auto layer = layers.begin(); layer != layers.end(); ++layer)
		{
			topIndex = BinarySearchInVectorSorted(mesh[2], layer->top);
			bottomIndex = BinarySearchInVectorSorted(mesh[2], layer->bottom);

			for (int iz = bottomIndex; iz < topIndex; iz++)
				for (int iy = 0; iy < mesh[1].size()-1; iy++)
				{
					cell = iz * nxy + iy * nx;
					for (int ix = 0; ix < nx; ix++, cell++)
						layerMaterials[cell] = layer->materialValue;
				}
		}
		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : SetLayersMaterials : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Apply objects materials to mesh
int SetObjectsMaterials(vector<GeoBox<double>> &objects, vector<Material> &objectMaterials, vector<double> *mesh)
{
	try
	{
		int ind[6], cell;
		int nx = mesh[0].size() - 1;
		int nxy = nx * (mesh[1].size() - 1);

		for (auto obj = objects.begin(); obj != objects.end(); ++obj)
		{
			for (int i = 0; i < 6; i++)
				ind[i] = BinarySearchInVectorSorted(mesh[i / 2], obj->coordinates[i]);

			for (int iz = ind[4]; iz < ind[5]; iz++)
				for (int iy = ind[2]; iy < ind[3]; iy++)
					for (int ix = ind[0]; ix < ind[1]; ix++)
						objectMaterials[iz * nxy + iy * nx + ix] = obj->materialValue;
		}
		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : SetObjectsMaterials : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Set primary field calculation area metarial parameters
int SetNormalMaterials(vector<Material> &layerMaterials, vector<Material> &objectMaterials)
{
	try
	{
		for (int materialIndex = 0; materialIndex < objectMaterials.size(); materialIndex++)
			if (!objectMaterials[materialIndex].SigmaH.Equal(layerMaterials[materialIndex].SigmaH))
				objectMaterials[materialIndex].SigmaN.constantValue = layerMaterials[materialIndex].SigmaH.constantValue;
		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : SetNormalMaterials : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Set material numbers to cells
int SetMaterialNumbers(vector<Material> &objectMaterials, vector<Material> &layerMaterials, vector<int> &materialNumbers, vector<Material> &uniqueMaterials, vector<double> &meshZ, vector<GeoLayer<double>> &layers, int nxy)
{
	vector<bool> done;
	done.assign(objectMaterials.size(), false);
	materialNumbers.resize(objectMaterials.size());
	try
	{
		int materialNumber = 1;

		// air
		int index1 = objectMaterials.size() - 1;
		materialNumbers[index1] = materialNumber;

		auto &material1 = objectMaterials[index1];
		uniqueMaterials.push_back(material1);

		for (int index2 = index1 - 1; index2 >= 0; index2--)
		{
			if (done[index2])
				continue;

			if (!material1.Equal(objectMaterials[index2]))
				continue;

			materialNumbers[index2] = materialNumber;
			done[index2] = true;
		}
		done[index1] = true;
		materialNumber++;

		// layers
		int zIndex = 0;
		int layerIndex = layers.size() - 1;
		for (index1 = 0; index1 < layerMaterials.size(); index1 += nxy, zIndex++)
		{
			if (done[index1])
				continue;

			auto &material1 = layerMaterials[index1];
			double z = (meshZ[zIndex] + (zIndex == meshZ.size()-1 ? meshZ[zIndex] : meshZ[zIndex+1])) * 0.5;
			while (z > layers[layerIndex].top)
			{
				materialNumber++;
				layerIndex--;

				uniqueMaterials.push_back(material1);
			}

			if (index1 == 0 && uniqueMaterials.size() == 1)
				uniqueMaterials.push_back(material1);

			for (int index2 = index1 + 1; index2 < index1 + nxy; index2++)
			{
				if (done[index2]) continue;
				if (!material1.Equal(objectMaterials[index2])) continue;

				materialNumbers[index2] = materialNumber;
				done[index2] = true;
			}

			if (material1.Equal(objectMaterials[index1]))
			{
				materialNumbers[index1] = materialNumber;
				done[index1] = true;
			}
		}

		materialNumber++;
		// other
		for (index1 = 0; index1 < objectMaterials.size(); index1++)
		{
			if (index1 == 7724)
				index1 = index1;
			if (done[index1])
				continue;

			materialNumbers[index1] = materialNumber;

			auto &material1 = objectMaterials[index1];
			uniqueMaterials.push_back(material1);

			for (int index2 = index1 + 1; index2 < objectMaterials.size(); index2++)
			{

				if (index2 == 7724)
					index2 = index2;

				if (done[index2])
					continue;

				if (!material1.Equal(objectMaterials[index2]))
					continue;

				materialNumbers[index2] = materialNumber;
				done[index2] = true;
			}
			done[index1] = true;
			materialNumber++;
		}
		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : SetMaterialNumbers : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Set materials to mesh
int SetMaterials(vector<GeoLayer<double>> &layers, vector<GeoBox<double>> &objects, vector<int> &materialNumbers, vector<double> *mesh, vector<Material> &uniqueMaterials)
{
	try
	{
		vector<Material> layerMaterials;
		vector<Material> objectMaterials;
		if (SetLayersMaterials(layers, layerMaterials, mesh) != 0) { write_to_log("Error : SetMaterials : Could not set layers materials to mesh\n"); return 1; }
		objectMaterials = layerMaterials;
		if (SetObjectsMaterials(objects, objectMaterials, mesh) != 0) { write_to_log("Error : SetMaterials : Could not set objects materials to mesh\n"); return 1; }
		if (SetNormalMaterials(layerMaterials, objectMaterials) != 0) { write_to_log("Error : SetMaterials : Could not set normal parameters to materials\n"); return 1; }
		if (SetMaterialNumbers(objectMaterials, layerMaterials, materialNumbers, uniqueMaterials, mesh[2], layers, (mesh[0].size() - 1) * (mesh[1].size() - 1)) != 0) { write_to_log("Error : SetMaterials : Could not set material numbers\n"); return 1; }

		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : SetMaterials : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Find edge by nodes
int FindEdge(int node1, int node2, map<pair<int, int>, int> &edgesToNumbers)
{
	static pair<int, int> edge;
	edge.first = node1;
	edge.second = node2;

	auto &found = edgesToNumbers.find(edge);
	if (found == edgesToNumbers.end())
	{
		char buf[2048];
		sprintf(buf, "Error : Could not find edge for nodes %d and %d\n", node1, node2);
		write_to_log(buf);
		exit(1);
	}
	else
		return found->second;
}

// Colelct mesh edges
int BuildEdges(vector<double> *mesh, set<pair<int, int>> &edges)
{
	try
	{
		pair<int, int> edge;
		int nx, ny, nz;
		nx = mesh[0].size();
		ny = mesh[1].size();
		nz = mesh[2].size();
		int nxy = nx*ny;
		nx = mesh[0].size();
		ny = mesh[1].size();
		nz = mesh[2].size();

		int cell = 0;
		int node;
		for (int iz = 0; iz < nz - 1; iz++)
		{
			for (int iy = 0; iy < ny - 1; iy++)
			{
				node = iz * nxy + iy * nx + 1;
				for (int ix = 0; ix < nx - 1; ix++, node++)
				{
					edges.insert(pair<int,int>(node, node + 1));
					edges.insert(pair<int,int>(node + nx, node + nx + 1));
					edges.insert(pair<int,int>(node + nxy, node + nxy + 1));
					edges.insert(pair<int,int>(node + nxy + nx, node + nxy + nx + 1));

					edges.insert(pair<int,int>(node, node + nx));
					edges.insert(pair<int,int>(node + 1, node + nx + 1));
					edges.insert(pair<int,int>(node + nxy, node + nxy + nx));
					edges.insert(pair<int,int>(node + nxy + 1, node + nxy + nx + 1));

					edges.insert(pair<int,int>(node, node + nxy));
					edges.insert(pair<int,int>(node + 1, node + 1 + nxy));
					edges.insert(pair<int,int>(node + nx, node + nx + nxy));
					edges.insert(pair<int,int>(node + nx + 1, node + nx + 1 + nxy));
				}
			}
		}
		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : SetObjectsMaterials : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Colelct mesh elements described with edges
int BuildElementsByEdges(vector<double> *mesh, set<pair<int, int>> &edges, vector<vector<int>> &elementsByEdges)
{
	try
	{
		pair<int, int> edge;
		int nx, ny, nz;
		nx = mesh[0].size();
		ny = mesh[1].size();
		nz = mesh[2].size();
		int nxy = nx*ny;
		nx = mesh[0].size();
		ny = mesh[1].size();
		nz = mesh[2].size();
		elementsByEdges.resize((nx - 1) * (ny - 1) * (nz - 1));
		map<pair<int, int>, int> edgesToNumbers;
		int number = 1;
		for (auto edge = edges.begin(); edge != edges.end(); ++edge, number++)
			edgesToNumbers[*edge] = number;

		int cell = 0;
		int node;
		for (int iz = 0; iz < nz - 1; iz++)
		{
			for (int iy = 0; iy < ny - 1; iy++)
			{
				node = iz * nxy + iy * nx + 1;
				for (int ix = 0; ix < nx - 1; ix++, node++, cell++)
				{
					elementsByEdges[cell].resize(12);

					elementsByEdges[cell][0] = FindEdge(node, node + 1, edgesToNumbers);
					elementsByEdges[cell][1] = FindEdge(node + nx, node + nx + 1, edgesToNumbers);
					elementsByEdges[cell][2] = FindEdge(node + nxy, node + nxy + 1, edgesToNumbers);
					elementsByEdges[cell][3] = FindEdge(node + nxy + nx, node + nxy + nx + 1, edgesToNumbers);

					elementsByEdges[cell][4] = FindEdge(node, node + nx, edgesToNumbers);					
					elementsByEdges[cell][5] = FindEdge(node + nxy, node + nxy + nx, edgesToNumbers);
					elementsByEdges[cell][6] = FindEdge(node + 1, node + nx + 1, edgesToNumbers);
					elementsByEdges[cell][7] = FindEdge(node + nxy + 1, node + nxy + nx + 1, edgesToNumbers);

					elementsByEdges[cell][8] = FindEdge(node, node + nxy, edgesToNumbers);
					elementsByEdges[cell][9] = FindEdge(node + 1, node + 1 + nxy, edgesToNumbers);
					elementsByEdges[cell][10] = FindEdge(node + nx, node + nx + nxy, edgesToNumbers);
					elementsByEdges[cell][11] = FindEdge(node + nx + 1, node + nx + 1 + nxy, edgesToNumbers);
				}
			}
		}
		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : SetObjectsMaterials : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Build mesh for 1D problem
int Build1DMesh(vector<double> &mesh1D)
{
	try
	{
		mesh1D.reserve(1000);

		double z = 0;
		double step = 1e-2;
		double sparse = 1.02;
		while (z > -1e+7)
		{
			mesh1D.push_back(z);
			z -= step;
			step *= sparse;
		}
		mesh1D.push_back(z);
		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : Build1DMesh : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Get total edges
int CalculateEdgesCount(vector<double> *mesh)
{
	try
	{
		int nx, ny, nz;
		int edgesCount = 0;
		nx = mesh[0].size();
		ny = mesh[1].size();
		nz = mesh[2].size();
		edgesCount += (nx - 1) * ny * nz;
		edgesCount += (ny - 1) * nx * nz;
		edgesCount += (nz - 1) * nx * ny;
		return edgesCount;
	}
	catch (exception ex)
	{
		write_to_log("Error : WriteTSize3DNode : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 0;
	}
}

// Mesh building procedure
int BuildMesh(vector<Receiver> &receivers, vector<GeoLayer<double>> &layers, vector<GeoBox<double>> &objects, Settings &meshSettings, vector<double> *meshes1D, vector<double> *meshesTemplate1D, vector<int> &materialNumbers, set<pair<int, int>> &edges, vector<vector<int>> &elementsByEdges, vector<Material> &uniqueMaterials)
{
	try
	{
		double *bounds[3] = { new double[4], new double[4], new double[4] };
		if (CalculateModelBounds(receivers, meshSettings, bounds) != 0) { write_to_log("Error : BuildMesh : Could not calculate model bounds\n"); return 1; }
		if (BuildBaseMesh3D(meshSettings, bounds, meshes1D) != 0) { write_to_log("Error : BuildMesh : Could not build base mesh\n"); return 1; }
		if (AlignMeshBounds(bounds, meshes1D) != 0) { write_to_log("Error : BuildMesh : Could not align base mesh\n"); return 1; }
		if (BuildBaseMesh3DTemplate(bounds, meshes1D, meshesTemplate1D) != 0) { write_to_log("Error : BuildMesh : Could not build template mesh\n"); return 1; }
		if (InsertObjectsIntoMesh(objects, meshSettings.eps, meshes1D) != 0) { write_to_log("Error : BuildMesh : Could not insert objects into mesh\n"); return 1; }
		if (InsertLayersIntoMesh(layers, meshSettings.eps, meshes1D) != 0) { write_to_log("Error : BuildMesh : Could not insert layers into mesh\n"); return 1; }
		if (SetMaterials(layers, objects, materialNumbers, meshes1D, uniqueMaterials) != 0) { write_to_log("Error : BuildMesh : Could not set material numbers to mesh\n"); return 1; }
		if (BuildEdges(meshes1D, edges) != 0) { write_to_log("Error : BuildMesh : Could not build edges\n"); return 1; }
		if (BuildElementsByEdges(meshes1D, edges, elementsByEdges) != 0) { write_to_log("Error : BuildMesh : Could not build elements by edges\n"); return 1; }

		delete[] bounds[0];
		delete[] bounds[1];
		delete[] bounds[2];
		return 0;
	}
	catch (exception ex)
	{
		write_to_log("Error : BuildMesh : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Read objects from file
int ReadObjects(char *file_name, vector<GeoBox<double>> &objects)
{
	FILE *file_in = NULL;
	try
	{
		if (!(file_in = fopen(file_name, "r"))) return 1;

		Material newMaterial;
		GeoBox<double> newObject;

		int objectsCount;
		fscanf(file_in, "%d", &objectsCount);
		objects.resize(objectsCount);

		for (int objectIndex = 0; objectIndex < objectsCount; objectIndex++)
		{
			fscanf(file_in, "%lf%lf%lf%lf%lf%lf%lf", newObject.coordinates, newObject.coordinates + 1, newObject.coordinates + 2, newObject.coordinates + 3, newObject.coordinates + 4, newObject.coordinates + 5, &newMaterial.SigmaH.constantValue);
			newMaterial.SigmaV.constantValue = newMaterial.SigmaH.constantValue;
			newMaterial.SigmaN.constantValue = newMaterial.SigmaH.constantValue;
			
			newMaterial.SigmaH.valueType = PhysicalValueType::Constant;
			newMaterial.SigmaV.valueType = PhysicalValueType::Constant;
			newMaterial.SigmaN.valueType = PhysicalValueType::Constant;
			newObject.materialValue = newMaterial;
			objects[objectIndex] = newObject;
		}

		fclose(file_in);
		return 0;
	}
	catch (exception ex)
	{
		if (file_in != NULL)
			fclose(file_in);
		write_to_log("Error : ReadObjects : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Create material for air
void CreateAirMaterial(Material &airMaterial)
{
	airMaterial.SigmaH.constantValue = 1e-8;
	airMaterial.SigmaV.constantValue = 1e-8;
	airMaterial.SigmaN.constantValue = 1e-8;
	airMaterial.SigmaH.valueType = PhysicalValueType::Constant;
	airMaterial.SigmaV.valueType = PhysicalValueType::Constant;
	airMaterial.SigmaN.valueType = PhysicalValueType::Constant;
}

// Read layers from file
int ReadZSig2D(char *file_name, vector<GeoLayer<double>> &layers)
{
	FILE *file_in = NULL;
	try
	{
		Material newMaterial;
		GeoLayer<double> newLayer;
		double z, sigma;
		int layersCount;
		if (!(file_in = fopen(file_name, "r"))) return 1;

		fscanf(file_in, "%d%", &layersCount);

		layers.resize(layersCount + 1);

		CreateAirMaterial(newMaterial);
		layers[0].materialValue = newMaterial;
		layers[0].bottom = 0;
		layers[0].top = 1e+30;
		layers[0].material = 0;
		for (int layerIndex = 0; layerIndex < layersCount; layerIndex++)
		{
			fscanf(file_in, "%lf%lf", &z, &sigma);
			newLayer.top = z;
			newLayer.bottom = layerIndex == 0 ? -1e+30 : z - layers[layers.size() - layerIndex - 2].top;
			newMaterial.SigmaH.constantValue = sigma;
			newMaterial.SigmaV.constantValue = sigma;
			newMaterial.SigmaN.constantValue = sigma;
			newMaterial.SigmaH.valueType = PhysicalValueType::Constant;
			newMaterial.SigmaV.valueType = PhysicalValueType::Constant;
			newMaterial.SigmaN.valueType = PhysicalValueType::Constant;
			newLayer.materialValue = newMaterial;
			layers[layers.size() - layerIndex - 1] = newLayer;
		}

		fclose(file_in);
		return 0;
	}
	catch (exception ex)
	{
		if (file_in != NULL)
			fclose(file_in);
		write_to_log("Error : ReadZSig2D : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Read receivers from file
int ReadReceivers(char *file_name, vector<Receiver> &receivers)
{
	FILE *file_in = NULL;
	try
	{
		if (!(file_in = fopen(file_name, "r"))) return 1;

		Receiver newReceiver;

		int receiversCount;
		fscanf(file_in, "%d", &receiversCount);
		receivers.resize(receiversCount);

		for (int receiverIndex = 0; receiverIndex < receiversCount; receiverIndex++)
		{
			fscanf(file_in, "%lf%lf%lf", newReceiver.coordinates, newReceiver.coordinates + 1, newReceiver.coordinates + 2);
			receivers[receiverIndex] = newReceiver;
		}

		fclose(file_in);
		return 0;
	}
	catch (exception ex)
	{
		if (file_in != NULL)
			fclose(file_in);
		write_to_log("Error : ReadReceivers : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Read settings from file
int ReadSettings(char *file_name, Settings &settings)
{
	FILE *file_in = NULL;
	try
	{
		if (!(file_in = fopen(file_name, "r"))) return 1;
		char c;

		fscanf(file_in, "%lf%lf", settings.steps, settings.steps + 2); settings.steps[1] = settings.steps[0];
		fscanf(file_in, "%lf%lf", settings.sparse, settings.sparse + 2); settings.sparse[1] = settings.sparse[0];
		fscanf(file_in, "%lf%lf", settings.gapFromReceivers, settings.gapFromReceivers + 2); settings.gapFromReceivers[1] = settings.gapFromReceivers[0];
		fscanf(file_in, "%lf%lf", settings.farBound, settings.farBound + 2); settings.farBound[1] = settings.farBound[0];
		fscanf(file_in, "%c", &c);
		while (c != '\n')
			fscanf(file_in, "%c", &c);
		fscanf(file_in, "%c", &c);
		if (c == 'x' || c == 'X')
			settings.currentDirection = 1;
		else
			if (c == 'y' || c == 'Y')
				settings.currentDirection = 0;
			else
			{
				write_to_log("Error : ReadSettings : Could not recognize current direction\n");
				fclose(file_in);
				return 1;
			}

		fscanf(file_in, "%lf", &settings.frequency);
		fscanf(file_in, "%d", &settings.threadsCount);
		fscanf(file_in, "%lf", &settings.SLAESolutionEps);
		fscanf(file_in, "%d", &settings.maxIterationsCount);

		fclose(file_in);
		return 0;
	}
	catch (exception ex)
	{
		if (file_in != NULL)
			fclose(file_in);
		write_to_log("Error : ReadSettings : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Write mesh file with nodes and cells count
int WriteInftry(char *file_name, int nx, int ny, int nz)
{
	FILE *file_out = NULL;
	try
	{
		if (!(file_out = fopen(file_name, "w"))) return 1;

		int nodesCount = nx * ny * nz;
		int cellsCount = (nx - 1) * (ny - 1) * (nz - 1);
		int bc1NodesCount = (nx * ny + nx * (nz - 2)  + (ny - 2) * (nz - 2)) * 2;

		fprintf(file_out, " ISLAU= %d", 0);
		fprintf(file_out, " INDKU1= %d", 0);
		fprintf(file_out, " INDFPO= %d", 0);
		fprintf(file_out, "\nKUZLOV= %d", nodesCount);
		fprintf(file_out, "   KPAR= %d", cellsCount);
		fprintf(file_out, "    KT1= %d", bc1NodesCount);
		fprintf(file_out, "   KTR2= %d", 0);
		fprintf(file_out, "   KTR3= %d\n", 0);
		fprintf(file_out, "   KT8= %d", 0);
		fprintf(file_out, "   KT9= %d\n", 0);
		fprintf(file_out, "KISRS1= %d", 0);
		fprintf(file_out, " KISRS2= %d", 0);
		fprintf(file_out, " KISRS3= %d", 0);
		fprintf(file_out, "   KBRS= %d\n", 0);
		fprintf(file_out, "   KT7= %d", 0);
		fprintf(file_out, "   KT10= %d", 0);
		fprintf(file_out, "  KTR4= %d", 0);
		fprintf(file_out, "  KTSIM= %d\n", 0);
		fprintf(file_out, "   KT6= %d", 0);
		fclose(file_out);

		return 0;
	}
	catch (exception ex)
	{
		if (file_out != NULL)
			fclose(file_out);
		write_to_log("Error : WriteInftry : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Write mesh nodes coordinates
int WriteXyz(char *file_name, vector<double> *mesh)
{
	FILE *file_out = NULL;
	try
	{
		double point[3];

		if (!(file_out = fopen(file_name, "wb"))) return 1;

		auto &meshX = mesh[0];
		auto &meshY = mesh[1];
		auto &meshZ = mesh[2];

		for (int iz = 0; iz < meshZ.size(); iz++)
		{
			point[2] = meshZ[iz];
			for (int iy = 0; iy < meshY.size(); iy++)
			{
				point[1] = meshY[iy];
				for (int ix = 0; ix < meshX.size(); ix++)
				{
					point[0] = meshX[ix];
					fwrite(point, sizeof(double), 3, file_out);
				}
			}
		}

		fclose(file_out);
		return 0;
	}
	catch (exception ex)
	{
		if (file_out != NULL)
			fclose(file_out);
		write_to_log("Error : WriteXyz : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Write mesh cells described with node numbers
int WriteNver(char *file_name, int nx, int ny, int nz)
{
	FILE *file_out = NULL;
	try
	{
		int element[14];
		element[8] = 0;
		element[9] = 0;
		element[10] = 0;
		element[11] = 0;
		element[12] = 0;
		element[13] = 0;
		int nxy = nx * ny;

		if (!(file_out = fopen(file_name, "wb"))) return 1;

		int node;
		for (int iz = 0; iz < nz - 1; iz++)
		{
			for (int iy = 0; iy < ny - 1; iy++)
			{
				node = nxy * iz + nx * iy + 1;
				element[0] = node;
				element[1] = node + 1;
				element[2] = node + nx;
				element[3] = node + nx + 1;
				element[4] = node + nxy;
				element[5] = node + nxy + 1;
				element[6] = node + nxy + nx;
				element[7] = node + nxy + nx + 1;
				for (int ix = 0; ix < nx - 1; ix++)
				{
					fwrite(element, sizeof(int), 14, file_out);

					element[0]++;
					element[1]++;
					element[2]++;
					element[3]++;
					element[4]++;
					element[5]++;
					element[6]++;
					element[7]++;
				}
			}
		}
		fclose(file_out);

		return 0;
	}
	catch (exception ex)
	{
		if (file_out != NULL)
			fclose(file_out);
		write_to_log("Error : WriteNver : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Write mesh cells material numbers
int WriteNvkat(char *file_name, vector<int> &materialNumbers)
{
	FILE *file_out = NULL;
	try
	{
		int m = 1;

		if (!(file_out = fopen(file_name, "wb"))) return 1;

		for (int i = 0; i < materialNumbers.size(); i++)
		{
			m = materialNumbers[i];
			fwrite(&m, sizeof(int), 1, file_out);
		}
		fclose(file_out);

		return 0;
	}
	catch (exception ex)
	{
		if (file_out != NULL)
			fclose(file_out);
		write_to_log("Error : WriteNvkat : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Get cell (with ix, iy, iz 1D mesh indeces) neighbours
void GetNeighbours(int nx, int ny, int nz, int ix, int iy, int iz, int cell, int *neighbours)
{
	int nxy = (nx-1) * (ny-1);
	neighbours[0] = ix == 0 ? -1 : 1 + cell - 1;
	neighbours[3] = ix == nx - 2 ? -1 : 1 + cell + 1;
	neighbours[1] = iy == 0 ? -1 : 1 + cell - (nx-1);
	neighbours[4] = iy == ny - 2 ? -1 : 1 + cell + (nx - 1);
	neighbours[2] = iz == 0 ? -1 : 1 + cell - nxy;
	neighbours[5] = iz == nz - 2 ? -1 : 1 + cell + nxy;
}

// Write each cell neighbours to file
int WriteElemNeib(char *file_name, int nx, int ny, int nz)
{
	FILE *file_out = NULL;
	try
	{
		int m = 0;
		int neighbours[6];

		if (!(file_out = fopen(file_name, "wb"))) return 1;

		int cell = 0;
		for (int iz = 0; iz < nz - 1; iz++)
		{
			for (int iy = 0; iy < ny - 1; iy++)
			{
				for (int ix = 0; ix < nx - 1; ix++, cell++)
				{
					GetNeighbours(nx, ny, nz, ix, iy, iz, cell, neighbours);
					for (int i = 0; i < 6; i++)
					{
						if (neighbours[i] == -1)
						{
							m = 0;
							fwrite(&m, sizeof(int), 1, file_out);
						}
						else
						{
							m = 1;
							fwrite(&m, sizeof(int), 1, file_out);
							fwrite(neighbours + i, sizeof(int), 1, file_out);
						}
					}
				}
			}
		}
		fclose(file_out);

		return 0;
	}
	catch (exception ex)
	{
		if (file_out != NULL)
			fclose(file_out);
		write_to_log("Error : WriteElemNeib : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Write material properties to file
int WriteProperty3D(char *fileName, vector<Material> &materials, PhysicalValueName prop1, PhysicalValueName prop2, bool revert = false)
{
	FILE *outputFile = NULL;
	try
	{
		int n;

		if (open_file_w(fileName, &outputFile) != 0)
			return 1;

		for (int i = 0; i < materials.size(); i++)
			if (!revert)
				fprintf(outputFile, "%d\t%.7e\t%.7e\n", i + 1, materials[i].GetValue(prop1).constantValue, materials[i].GetValue(prop2).constantValue);
			else
				fprintf(outputFile, "%d\t%.7e\t%.7e\n", i + 1, 1.0 / materials[i].GetValue(prop1).constantValue, 1.0 / materials[i].GetValue(prop2).constantValue);

		fclose(outputFile);
		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		write_to_log("Error : WriteProperty3D : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Write empty nodal T-Matrix size
int WriteTSize3DNode(char *fileName, vector<double> *mesh)
{
	FILE *outputFile = NULL;
	try
	{
		try
		{
			if (open_file_w(fileName, &outputFile) != 0)
				return 1;
			
			fprintf(outputFile, "%d\n", 0);

			fclose(outputFile);
			return 0;
		}
		catch (exception ex)
		{
			if (outputFile != NULL)
				fclose(outputFile);
			return 1;
		}
		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		write_to_log("Error : WriteTSize3DNode : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Write empty edge T-Matrix size
int WriteTSize3DEdge(char *fileName, int edgesCount)
{
	FILE *outputFile = NULL;
	try
	{
		if (open_file_w(fileName, &outputFile) != 0)
			return 1;

		fprintf(outputFile, "%d\n%d\n", 0, edgesCount);

		fclose(outputFile);
		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		write_to_log("Error : WriteTSize3DNode : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Write 1D meshes to file
int Write3DMeshRegular(char *fileName, vector<double> *mesh)
{
	FILE *outputFile = NULL;
	try
	{
		if (open_file_w(fileName, &outputFile) != 0)
			return 1;

		int nx, ny, nz;
		int edgesCount = 0;
		nx = mesh[0].size();
		ny = mesh[1].size();
		nz = mesh[2].size();
		edgesCount += (nx - 1) * ny * nz;
		edgesCount += (ny - 1) * nx * nz;
		edgesCount += (nz - 1) * nx * ny;
		for (int direction = 0; direction < 3; direction++)
		{
			fprintf(outputFile, "%d\n", mesh[direction].size());
			for (int meshIndex = 0; meshIndex < mesh[direction].size(); meshIndex++)
			{
				fprintf(outputFile, "%d\t%.14e\n", meshIndex + 1, mesh[direction][meshIndex]);
			}
		}

		fprintf(outputFile, "0\n0\n");

		fclose(outputFile);
		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		write_to_log("Error : Write3DMeshRegular : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Write regular mesh cells flags
int WriteRegular(char *fileName, int nx, int ny, int nz)
{
	FILE *outputFile = NULL;
	try
	{
		int m = 0;
		int neighbours[6];

		if (!(outputFile = fopen(fileName, "wb"))) return 1;

		int cellsCount = (nx - 1) * (ny - 1) * (nz - 1);
		fwrite(&cellsCount, sizeof(int), 1, outputFile);
		for (int cell = 1; cell <= cellsCount; cell++)
			fwrite(&cell, sizeof(int), 1, outputFile);
		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		write_to_log("Error : WriteRegular : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Write boundary nodes
int WriteL13D(char *fileName, int nx, int ny, int nz)
{
	FILE *outputFile = NULL;
	try
	{
		int m = 0;
		int neighbours[6];

		if (!(outputFile = fopen(fileName, "wb"))) return 1;

		int node = 1;
		for (int iz = 0; iz < nz; iz++)
			for (int iy = 0; iy < ny; iy++)
				for (int ix = 0; ix < nx; ix++, node++)
					if (ix == 0 || ix == nx - 1 || iy == 0 || iy == ny - 1 || iz == 0 || iz == nz - 1)
						fwrite(&node, sizeof(int), 1, outputFile);
		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		write_to_log("Error : WriteL13D : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Write node numbers for each edge
int WriteNodesForEdges(char *fileName, set<pair<int, int>> &edges)
{
	FILE *outputFile = NULL;
	try
	{
		if (!(outputFile = fopen(fileName, "wb"))) return 1;

		for (auto edge = edges.begin(); edge != edges.end(); ++edge)
		{
			fwrite(&edge->first, sizeof(int), 1, outputFile);
			fwrite(&edge->second, sizeof(int), 1, outputFile);
		}
		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		write_to_log("Error : WriteNodesForEdges : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Write mesh cells described by edges
int WriteEdges(char *fileName, vector<vector<int>> &elementsByEdges)
{
	FILE *outputFile = NULL;
	try
	{
		int zero = 0;
		int one = 1;
		if (!(outputFile = fopen(fileName, "wb"))) return 1;

		for (auto element = elementsByEdges.begin(); element != elementsByEdges.end(); ++element)
		{
			for (int i = 0; i < element->size(); i++)
				fwrite(&((*element)[i]), sizeof(int), 1, outputFile);

			for (int i = 0; i < 12; i++)
				fwrite(&zero, sizeof(int), 1, outputFile);
			fwrite(&one, sizeof(int), 1, outputFile);
		}
		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		write_to_log("Error : WriteNodesForEdges : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Write empty edge T-Matrix
int WriteIg3D(char *fileName)
{
	FILE *outputFile = NULL;
	try
	{		
		if (!(outputFile = fopen(fileName, "wb"))) return 1;
		int m = 1;
		fwrite(&m, sizeof(int), 1, outputFile);
		fwrite(&m, sizeof(int), 1, outputFile);
		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		write_to_log("Error : WriteIg3D : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}
int WriteJg3D(char *fileName)
{
	FILE *outputFile = NULL;
	try
	{
		if (!(outputFile = fopen(fileName, "wb"))) return 1;
		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		write_to_log("Error : WriteJg3D : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}
int WriteGg3D(char *fileName)
{
	FILE *outputFile = NULL;
	try
	{
		if (!(outputFile = fopen(fileName, "wb"))) return 1;
		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		write_to_log("Error : WriteGg3D : ");
		write_to_log(ex.what());
		write_to_log("\n");
		return 1;
	}
}

// Write current direction file
int WriteAlpha(char *fileName, Settings &meshSettings)
{
	FILE *outputFile = NULL;
	try
	{
		if (!(outputFile = fopen(fileName, "w"))) return 1;
		fprintf(outputFile, "%d\n", meshSettings.currentDirection);
		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		return 1;
	}
}

// Write frequency file
int WriteNu(char *fileName, Settings &meshSettings)
{
	FILE *outputFile = NULL;
	try
	{
		if (!(outputFile = fopen(fileName, "w"))) return 1;
		fprintf(outputFile, "%.13e\n%d", meshSettings.frequency, 3);
		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		return 1;
	}
}

// Write receivers coordinates
int WritePointRes(char *fileName, vector<Receiver> &receivers)
{
	FILE *outputFile = NULL;
	try
	{
		if (!(outputFile = fopen(fileName, "w"))) return 1;

		fprintf(outputFile, "%d\n", receivers.size());
		for (auto rec = receivers.begin(); rec != receivers.end(); ++rec)
			fprintf(outputFile, "%.13e\t%.13e\t%.13e\n", rec->coordinates[0], rec->coordinates[1], rec->coordinates[2]);

		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		return 1;
	}
}

// Write receivers coordinates
int WritexyzVectorE(char *fileName, vector<Receiver> &receivers)
{
	FILE *outputFile = NULL;
	try
	{
		if (!(outputFile = fopen(fileName, "w"))) return 1;

		fprintf(outputFile, "%d\n", receivers.size());
		for (auto rec = receivers.begin(); rec != receivers.end(); ++rec)
			fprintf(outputFile, "%.13e\t%.13e\t%.13e\n", rec->coordinates[0], rec->coordinates[1], rec->coordinates[2]);

		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		return 1;
	}
}

// Write 1D mesh file
int WriteMesh1D(char *fileName, vector<double> &mesh1D)
{
	FILE *outputFile = NULL;
	try
	{
		if (!(outputFile = fopen(fileName, "w"))) return 1;

		fprintf(outputFile, "\n\n");
		for (auto coord = mesh1D.rbegin(); coord != mesh1D.rend(); ++coord)
			fprintf(outputFile, "%.13e\n", *coord);

		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		return 1;
	}
}

// Write 1D problem mesh size file
int WriteSetka1DEy(char *fileName, vector<double> &mesh1D)
{
	FILE *outputFile = NULL;
	try
	{
		if (!(outputFile = fopen(fileName, "w"))) return 1;
		fprintf(outputFile, "%d\n", mesh1D.size());
		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		return 1;
	}
}

// Write 1D problem environment file
int WriteSreda1DAy(char *fileName, vector<GeoLayer<double>> &layers)
{
	FILE *outputFile = NULL;
	try
	{
		if (!(outputFile = fopen(fileName, "w"))) return 1;
		fprintf(outputFile, "%d\n", layers.size()-1);
		for (int layerIndex = layers.size() - 1; layerIndex > 0; layerIndex--)
			fprintf(outputFile, "%.13e\t%.13e\t%d\n", layers[layerIndex].top, layers[layerIndex].materialValue.SigmaH.constantValue, 1);

		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		return 1;
	}
}

// Write 3D problem SLAE size
int WriteKuSlau2(char *fileName, int edgesCount)
{
	FILE *outputFile = NULL;
	try
	{
		if (!(outputFile = fopen(fileName, "w"))) return 1;
		fprintf(outputFile, "%d\t%d\t%d\n", edgesCount * 2, 0, 0);
		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		return 1;
	}
}

// Write threads count file
int WriteNThreads(char *fileName, int threadsCount)
{
	FILE *outputFile = NULL;
	try
	{
		if (!(outputFile = fopen(fileName, "w"))) return 1;
		fprintf(outputFile, "%d\n", threadsCount);
		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		return 1;
	}
}

// Write SLAE solving parameters file
int WriteISSCTxt(char *fileName, double SLAEEps, int SLAEMaxIter)
{
	FILE *outputFile = NULL;
	try
	{
		if (!(outputFile = fopen(fileName, "w"))) return 1;

		fprintf(outputFile, "%d\n", SLAEMaxIter);
		fprintf(outputFile, "%.13e\n", SLAEEps);

		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		return 1;
	}
}

// Write simple needed file
int WriteHarm3DParams(char *fileName)
{
	FILE *outputFile = NULL;
	try
	{
		if (!(outputFile = fopen(fileName, "w"))) return 1;

		fprintf(outputFile, "0\n0\n");
		fclose(outputFile);

		return 0;
	}
	catch (exception ex)
	{
		if (outputFile != NULL)
			fclose(outputFile);
		return 1;
	}
}

// Write all mesh files procedure
int WriteMesh(vector<double> *meshes1D, vector<double> *meshesTemplate1D, vector<int> &materialNumbers, set<pair<int, int>> &edges, vector<vector<int>> &elementsByEdges, Settings &meshSettings)
{
	vector<double> mesh1D;
	int nx, ny, nz;
	int edgesCount;

	nx = meshes1D[0].size();
	ny = meshes1D[1].size();
	nz = meshes1D[2].size();
	edgesCount = CalculateEdgesCount(meshes1D);

	if (WriteInftry("inftry.dat", nx, ny, nz) != 0) return 1;
	if (WriteXyz("xyz.dat", meshes1D) != 0) return 1;
	if (WriteNver("nver.dat", nx, ny, nz) != 0) return 1;
	if (WriteNvkat("nvkat.dat", materialNumbers) != 0) return 1;
	if (WriteElemNeib("elem_neib", nx, ny, nz) != 0) return 1;
	if (WriteTSize3DNode("tsize3d.dat", meshes1D) != 0) return 1;
	if (WriteTSize3DEdge("tsize3d_.dat", edgesCount) != 0) return 1;
	if (Write3DMeshRegular("3dmeshregular", meshes1D) != 0) return 1;
	if (WriteRegular("regular", nx, ny, nz) != 0) return 1;
	if (WriteL13D("L13d.dat", nx, ny, nz) != 0) return 1;
	if (WriteNodesForEdges("nodesforedges.dat", edges) != 0) return 1;
	if (WriteEdges("edges.dat", elementsByEdges) != 0) return 1;
	if (WriteIg3D("ig3d_.dat") != 0) return 1;
	if (WriteJg3D("jg3d_.dat") != 0) return 1;
	if (WriteGg3D("gg3d_.dat") != 0) return 1;
	if (WriteAlpha("alfa", meshSettings) != 0) return 1;
	if (WriteNu("nu", meshSettings) != 0) return 1;
	if (Build1DMesh(mesh1D) != 0) return 1;
	if (WriteMesh1D("mesh1d", mesh1D) != 0) return 1;
	if (WriteSetka1DEy("setka1DEy", mesh1D) != 0) return 1;
	if (WriteKuSlau2("kuslau2", edgesCount) != 0) return 1;
	if (WriteNThreads("nthreads.txt", meshSettings.threadsCount) != 0) return 1;
	if (WriteISSCTxt("issc.txt", meshSettings.SLAESolutionEps, meshSettings.maxIterationsCount) != 0) return 1;
	if (WriteHarm3DParams("harm3dparams") != 0) return 1;

	return 0;
}

// Write material parameters files
int WriteMaterials(vector<Material> &materials)
{
	if (WriteProperty3D("Sig3d", materials, SigmaH, SigmaN, false) != 0)
	{
		write_to_log("Could not write 'Sig3d'\n");
		return 1;
	}

	if (WriteProperty3D("dpr3D", materials, Eps, Eps) != 0)
	{
		write_to_log("Could not write 'dpr3D'\n");
		return 1;
	}

	if (WriteProperty3D("mu3D", materials, Mu, Mu) != 0)
	{
		write_to_log("Could not write 'mu3D'\n");
		return 1;
	}
}

// Main mesh building and files writing procedure
int MainProcedure()
{
	char buf[2048];
	vector<Receiver> receivers;
	vector<GeoLayer<double>> layers;
	vector<GeoBox<double>> objects;
	vector<Material> materials;
	Settings meshSettings;
	vector<double> meshes1D[3];
	vector<double> meshesTemplate1D[3];
	vector<int> materialNumbers;
	set<pair<int, int>> edges;
	vector<vector<int>> elementsByEdges;

	fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
	write_to_log("Reading settings\n");
	if (ReadSettings("RegularMeshBuilderSettings.cfg", meshSettings) != 0)
	{
		write_to_log("Error : MainProcedure : Could not read RegularMeshBuilderSettings\n");
		return 1;
	}

	fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
	write_to_log("Reading layers\n");

	if (ReadZSig2D("z_sig_2d", layers) != 0)
	{
		write_to_log("Error : MainProcedure : Could not read z_sig_2d\n");
		return 1;
	}

	fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
	write_to_log("Reading objects\n");
	if (ReadObjects("objects", objects) != 0)
	{
		write_to_log("Error : MainProcedure : Could not read objects\n");
		return 1;
	}

	fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
	write_to_log("Reading receivers\n");
	if (ReadReceivers("xyzVectorB", receivers) != 0)
	{
		write_to_log("Error : MainProcedure : Could not read xyzVectorB\n");
		return 1;
	}

	fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
	write_to_log("Building mesh\n");
	if (BuildMesh(receivers, layers, objects, meshSettings, meshes1D, meshesTemplate1D, materialNumbers, edges, elementsByEdges, materials) != 0)
	{
		write_to_log("Error : MainProcedure : Could not build mesh\n");
		return 1;
	}

	fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
	write_to_log("Writing materials\n");
	if (WriteMaterials(materials) != 0)
	{
		write_to_log("Error : MainProcedure : Could not write materials\n");
		return 1;
	}
	
	fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
	write_to_log("Writing receivers\n");
	if (WritePointRes("pointres", receivers) != 0)
	{
		write_to_log("Error : MainProcedure : Could not write pointres\n");
		return 1;
	}

	fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
	if (WritexyzVectorE("xyzVectorE", receivers) != 0)
	{
		write_to_log("Error : MainProcedure : Could not write xyzVectorE\n");
		return 1;
	}

	
	fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
	write_to_log("Writing 1d layers\n");
	if (WriteSreda1DAy("sreda1D.ay", layers) != 0)
	{
		write_to_log("Error : MainProcedure : Could not write sreda1D.ay\n");
		return 1;
	}

	fprintf(log_file, "%lf :", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
	write_to_log("Writing mesh\n");
	if (WriteMesh(meshes1D, meshesTemplate1D, materialNumbers, edges, elementsByEdges, meshSettings) != 0)
	{
		write_to_log("Error : MainProcedure : Could not write mesh\n");
		return 1;
	}

	return 0;
}

int main(int argc, char **argv)
{
	tProgramStart = std::chrono::high_resolution_clock::now();
	char buf[2048];

	if (open_log("RegularMeshBuilder.log") != 0)
	{
		printf("Cound not open RegularMeshBuilder.log\n");
		return 1;
	}
	
	int status = MainProcedure();
	sprintf(buf, "Status = %d\n", status);
	write_to_log(buf);

	fprintf(log_file, "%lf : All done\n", std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - tProgramStart).count());
	fclose(log_file);
	return status;
}