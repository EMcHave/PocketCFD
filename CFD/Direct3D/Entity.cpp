#include "pch.h"
#include "Entity.h"

Entity::Entity()
{

}

Entity::Entity(
	DirectX::XMFLOAT4 size,
	InitialConditions initialConditions,
	Material material,
	EntityType type)
{
	m_size = size;
	m_position = initialConditions.position;
	m_material = material;

	switch (type)
	{
	case PlaneTriangle:
		CreateTriangularMesh1D();
		break;
	case PlaneQuad:
		break;
	case SolidQuad:
		CreateOrthogonalMesh();
		break;
	case Cylinder:
		CreateCylinderMesh();
		break;
	case Sphere:
		break;
	default:
		break;
	}
		
	SetInitialConditions(initialConditions.velocity);
}

void Entity::SetInitialConditions(DirectX::XMFLOAT3 velocity)
{
	for (auto particle : m_particles)
		particle->Velocity(velocity);
}

void Entity::CreateOrthogonalMesh()
{
	XMFLOAT3 n1(1, 0, 0);
	XMFLOAT3 n2(0, 1, 0);
	XMFLOAT3 n3(0, 0, 1);

	for (int k = 0; k < m_size.z; k++)
		for (int j = 0; j < m_size.y; j++)
			for (int i = 0; i < m_size.x; i++)
			{
				auto cube = std::make_shared<Cube>(XMFLOAT3(i * m_size.w + m_position.x, k * m_size.w + m_position.y, -j * m_size.w + m_position.z), m_size.w / 2,  false);
				m_particles.push_back(cube);
			}

	for (int k = 0; k < m_size.z; k++)
	{
		for (int j = 0; j < m_size.y; j++)
			for (int i = 0; i < m_size.x - 1; i++)
			{
				std::array<std::shared_ptr<Line>, 3> vectors1
				{
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n1)),
						std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n3)),
						std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n2))
				};
				std::array<std::shared_ptr<Line>, 3> vectors2
				{
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), -XMLoadFloat3(&n1)),
						std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n3)),
						std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n2))
				};
				auto connection = std::make_shared<Connection>(
					m_particles[k * m_size.x * m_size.y + j * m_size.x + i], m_particles[k * m_size.x * m_size.y + j * m_size.x + i + 1],
					vectors1, vectors2
				);

				m_connections.push_back(connection);
			}

		for (int j = 0; j < m_size.y - 1; j++)
			for (int i = 0; i < m_size.x; i++)
			{
				std::array<std::shared_ptr<Line>, 3> vectors1
				{
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), -XMLoadFloat3(&n3)),
						std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n1)),
						std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n2))
				};
				std::array<std::shared_ptr<Line>, 3> vectors2
				{
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n3)),
						std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n1)),
						std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n2))
				};
				auto connection = std::make_shared<Connection>(
					m_particles[k * m_size.x * m_size.y + j * m_size.x + i], m_particles[k * m_size.x * m_size.y + (j + 1) * m_size.x + i],
					vectors1, vectors2
				);
				m_connections.push_back(connection);
			}
	}

	for (int k = 0; k < m_size.z - 1; k++)
		for (int j = 0; j < m_size.y; j++)
			for (int i = 0; i < m_size.x; i++)
			{
				std::array<std::shared_ptr<Line>, 3> vectors1
				{
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n2)),
						std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n3)),
						std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), -XMLoadFloat3(&n1))
				};
				std::array<std::shared_ptr<Line>, 3> vectors2
				{
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), -XMLoadFloat3(&n2)),
						std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n3)),
						std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), -XMLoadFloat3(&n1))
				};
				auto connection = std::make_shared<Connection>(
					m_particles[k * m_size.x * m_size.y + j * m_size.x + i], m_particles[(k + 1) * m_size.x * m_size.y + j * m_size.x + i],
					vectors1, vectors2
				);
				m_connections.push_back(connection);
			}
}

void Entity::CreateTriangularMesh1D()
{
	XMFLOAT3 n1(1, 0, 0);
	XMFLOAT3 n2(0, 1, 0);
	XMFLOAT3 n3(0, 0, 1);

	XMVECTOR quat60 = XMQuaternionRotationNormal(XMLoadFloat3(&n2), XM_PI / 3);
	XMVECTOR quat120 = XMQuaternionRotationNormal(XMLoadFloat3(&n2), 2 * XM_PI / 3);

	XMVECTOR rotated = XMVector3Rotate(XMLoadFloat3(&n1), quat120);


	for (int j = 0; j < m_size.y; j++)
		for (int i = 0; i < m_size.x; i++)
		{
			std::shared_ptr<Cube> cube;
			if (j % 2 == 0)
			{
				cube = std::make_shared<Cube>(XMFLOAT3(
					i * m_size.w + m_position.x,
					m_position.y,
					-j * m_size.w * sinf(XM_PI / 3) + m_position.z), m_size.w * 0.75,  false);
			}
			else
			{
				cube = std::make_shared<Cube>(XMFLOAT3(
					i * m_size.w + m_size.w * cosf(XM_PI / 3) + m_position.x,
					m_position.y,
					-j * m_size.w * sinf(XM_PI / 3) + m_position.z), m_size.w * 0.75, false);
			}
				
			m_particles.push_back(cube);
		}

	for (int j = 0; j < m_size.y; j++)
		for (int i = 0; i < m_size.x - 1; i++)
		{
			std::array<std::shared_ptr<Line>, 3> vectors1
			{
				std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n1)),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n3)),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n2))
			};
			std::array<std::shared_ptr<Line>, 3> vectors2
			{
				std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), -XMLoadFloat3(&n1)),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n3)),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), XMLoadFloat3(&n2))
			};
			auto connection = std::make_shared<Connection>(
				m_particles[j * m_size.x + i], m_particles[j * m_size.x + i + 1],
				vectors1, vectors2
			);

			m_connections.push_back(connection);
		}

	for (int j = 0; j < m_size.y - 1; j++)
		for (int i = 0; i < m_size.x; i++)
		{
			auto v1 = (j % 2 == 0) * XMVector3Rotate(XMLoadFloat3(&n1), quat60) + (j % 2 != 0) * XMVector3Rotate(XMLoadFloat3(&n1), quat120);
			auto v2 = (j % 2 == 0) * XMVector3Rotate(XMLoadFloat3(&n3), quat60) + (j % 2 != 0) * XMVector3Rotate(XMLoadFloat3(&n3), quat120);
			auto v3 = XMLoadFloat3(&n3);

			std::array<std::shared_ptr<Line>, 3> vectors1
			{
				std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v1),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v2),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v3)
			};
			std::array<std::shared_ptr<Line>, 3> vectors2
			{
				std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), -v1),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v2),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v3)
			};
			auto connection = std::make_shared<Connection>(
				m_particles[j * m_size.x + i], m_particles[(j + 1) * m_size.x + i],
				vectors1, vectors2
			);
			m_connections.push_back(connection);
		}

	for (int i = 0; i < m_size.x - 1; i++)
		for (int j = 0; j < m_size.y - 1; j++)
		{
			auto v1_b = XMVector3Rotate(XMLoadFloat3(&n1), quat120);
			auto v2_b = XMVector3Rotate(XMLoadFloat3(&n3), quat120);

			auto v1_f = XMVector3Rotate(XMLoadFloat3(&n1), quat60);
			auto v2_f = XMVector3Rotate(XMLoadFloat3(&n3), quat60);

			auto v3 = XMLoadFloat3(&n3);

			std::array<std::shared_ptr<Line>, 3> vectors1_b
			{
				std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v1_b),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v2_b),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v3)
			};
			std::array<std::shared_ptr<Line>, 3> vectors2_b
			{
				std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), -v1_b),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v2_b),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v3)
			};

			std::array<std::shared_ptr<Line>, 3> vectors1_f
			{
				std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v1_f),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v2_f),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v3)
			};
			std::array<std::shared_ptr<Line>, 3> vectors2_f
			{
				std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), -v1_f),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v2_f),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v3)
			};

			if (j % 2 == 0)
			{
				auto connection = std::make_shared<Connection>(
					m_particles[j * m_size.x + (i + 1)], m_particles[(j + 1) * m_size.x + i],
					vectors1_b, vectors2_b
				);
				m_connections.push_back(connection);
			}
			else
			{
				auto connection2 = std::make_shared<Connection>(
					m_particles[j * m_size.x + i], m_particles[(j + 1) * m_size.x + (i + 1)],
					vectors1_f, vectors2_f
				);
				m_connections.push_back(connection2);
			}
		}

	for (auto particle : m_particles)
		if (particle->NumberOfConnections() < 6)
		{
			m_contactParticles.push_back(particle);
			particle->MakeContact();
		}
}

void Entity::CreateCylinderMesh()
{
	float r = m_size.x;
	float dfi = m_size.w / r;
	float NFI = 2 * XM_PI / dfi;

	XMFLOAT3 n1(1, 0, 0);
	XMFLOAT3 n2(0, 1, 0);
	XMFLOAT3 n3(0, 0, 1);

	for (int k = 0; k < m_size.z; k++)
		for (int j = 0; j < NFI; j++)
		{
			auto cube = std::make_shared<Cube>(XMFLOAT3(
				r * cosf(dfi * j) + m_position.x,
				k * m_size.w + m_position.y,
				r * sinf(dfi * j) + m_position.z), m_size.w * 0.75, true);
			XMVECTOR quat = XMQuaternionRotationNormal(XMLoadFloat3(&n2), -j * dfi);
			cube->SetInitialQuaternion(quat);
			m_particles.push_back(cube);
		}

	for (auto p : m_particles)
	{
		if (abs(p->Position().y) < m_size.w * 2)
		{
			m_contactParticles.push_back(p);
			p->MakeContact();
		}
	}
		

	for (int k = 0; k < m_size.z; k++)
	{
		for (int j = 0; j < NFI; j++)
		{
			std::shared_ptr<Cube> p2;
			auto p1 = m_particles[k * (int)(NFI + 1) + j];
			if (j + 1 == (int)NFI + 1)
				p2 = m_particles[k * (int)(NFI + 1)];
			else
				p2 = m_particles[k * (int)(NFI + 1) + j + 1];


			auto v1 = (p2->VectorPosition() - p1->VectorPosition()) / XMVector3Length(p2->VectorPosition() - p1->VectorPosition());
			auto v3 = XMLoadFloat3(&n2);
			auto v2 = XMVector3Cross(v1, v3);

			std::array<std::shared_ptr<Line>, 3> vectors1
			{
				std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v1),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v2),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v3)
			};
			std::array<std::shared_ptr<Line>, 3> vectors2
			{
				std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), -v1),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v2),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v3)
			};

			auto connection = std::make_shared<Connection>(
				p1, p2,
				vectors1, vectors2
			);

			m_connections.push_back(connection);
		}
	}

	for (int j = 0; j < NFI; j++)
		for (int k = 0; k < m_size.z - 1; k++)
		{
			auto p1 = m_particles[k * (int)(NFI + 1) + j];
			auto p2 = m_particles[(k + 1) * (int)(NFI + 1) + j];

			XMFLOAT3 z(0, k * m_size.w, 0);

			auto v1 = XMLoadFloat3(&n2);
			auto v2 = (XMLoadFloat3(&z) - p1->VectorPosition()) / XMVector3Length(XMLoadFloat3(&z) - p1->VectorPosition());
			auto v3 = XMVector3Cross(v2, v1);

			std::array<std::shared_ptr<Line>, 3> vectors1
			{
				std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v1),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v2),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v3)
			};
			std::array<std::shared_ptr<Line>, 3> vectors2
			{
				std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), -v1),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v2),
					std::make_shared<Line>(Axis::W, XMFLOAT3(0, 0, 0), v3)
			};

			auto connection = std::make_shared<Connection>(
				p1, p2,
				vectors1, vectors2
			);

			m_connections.push_back(connection);
		}
}

void Entity::RigidRotate(XMVECTOR axis, float angle_rad)
{
	auto quat = XMQuaternionRotationNormal(axis / XMVector3Length(axis), angle_rad);
	for (auto p : m_particles)
	{
		p->Quaternion(quat);
		p->Position(XMVector3Rotate(p->VectorPosition(), quat));
	}
	for (auto c : m_connections)
		c->Update();
}
