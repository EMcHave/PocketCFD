#pragma once


#include <ppltasks.h>
#include "Common/DeviceResources.h"
 //#include "Common/BasicReaderWriter.h"
#include "Common/BasicLoader.h"
#include "ConstantBuffers.h"
#include "LogicClass.h"
#include "Line.h"


using namespace winrt::Windows::Foundation::Numerics;



class MyRenderer : public std::enable_shared_from_this<MyRenderer>
{
public:
    MyRenderer(std::shared_ptr<DX::DeviceResources> const& deviceResources);

    void CreateDeviceDependentResources(_In_ std::shared_ptr<LogicClass> logic);
    void CreateWindowSizeDependentResources();
    void ReleaseDeviceDependentResources();
    void Render();


private:

    // Cached pointer to device resources.
    std::shared_ptr<DX::DeviceResources>            m_deviceResources;
    std::vector<std::shared_ptr<Cube>>              m_objects;
    std::vector<std::shared_ptr<Line>>              m_lines;
    std::shared_ptr <LogicClass>                    m_game;

    // Constant Buffers
    winrt::com_ptr<ID3D11VertexShader>          m_vertexShader;
    winrt::com_ptr<ID3D11PixelShader>           m_pixelShader;
    winrt::com_ptr<ID3D11InputLayout>           m_inputLayout;

    winrt::com_ptr<ID3D11Buffer>                m_vertexBuffer;
    winrt::com_ptr<ID3D11Buffer>                m_indexBuffer;
    winrt::com_ptr<ID3D11Buffer>                m_constantBufferChangeOnResize;
    winrt::com_ptr<ID3D11Buffer>                m_constantBufferChangesEveryFrame;
    winrt::com_ptr<ID3D11Buffer>                m_constantBufferChangesEveryPrim;


    bool	m_loadingComplete;
    float	m_degreesPerSecond;
    float   degree;
    bool	m_tracking;

    DirectX::XMMATRIX view;
    DirectX::XMFLOAT4X4 projection;
};

