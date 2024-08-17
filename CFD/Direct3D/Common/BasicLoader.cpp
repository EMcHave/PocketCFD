// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
//
// Copyright (c) Microsoft Corporation. All rights reserved

#include "pch.h"
#include "BasicLoader.h"

BasicLoader::BasicLoader(
    _In_ winrt::com_ptr<ID3D11Device3> d3dDevice,
    _In_opt_ IWICImagingFactory2* wicFactory
    ) :
    m_d3dDevice(std::move(d3dDevice))
{
    m_wicFactory.copy_from(wicFactory);
}

std::wstring BasicLoader::GetExtension(
    _In_ winrt::hstring const& filename
    )
{
    std::wstring extension;
    auto lastDot = std::find(filename.rbegin(), filename.rend(), L'.');
    if (lastDot != filename.rend())
    {
        std::transform(lastDot.base(), filename.end(), std::back_inserter(extension),
            [](wchar_t c) { return std::towlower(c); });
    }
    return extension;
}



void BasicLoader::CreateInputLayout(
    _In_reads_bytes_(bytecodeSize) const byte* bytecode,
    _In_ uint32_t bytecodeSize,
    _In_reads_opt_(layoutDescNumElements) const D3D11_INPUT_ELEMENT_DESC* layoutDesc,
    _In_ uint32_t layoutDescNumElements,
    _Out_ ID3D11InputLayout** layout
    )
{
    if (layoutDesc == nullptr)
    {
        // If no input layout is specified, use the BasicVertex layout.
        static constexpr D3D11_INPUT_ELEMENT_DESC basicVertexLayoutDesc[] =
        {
            { "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0,  D3D11_INPUT_PER_VERTEX_DATA, 0 },
            { "NORMAL",   0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0 },
            { "TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT,    0, 24, D3D11_INPUT_PER_VERTEX_DATA, 0 },
        };

        winrt::check_hresult(
            m_d3dDevice->CreateInputLayout(
                basicVertexLayoutDesc,
                ARRAYSIZE(basicVertexLayoutDesc),
                bytecode,
                bytecodeSize,
                layout
                )
            );
    }
    else
    {
        winrt::check_hresult(
            m_d3dDevice->CreateInputLayout(
                layoutDesc,
                layoutDescNumElements,
                bytecode,
                bytecodeSize,
                layout
                )
            );
    }
}





void BasicLoader::LoadShader(
    _In_ winrt::hstring const& filename,
    _In_reads_opt_(layoutDescNumElements) D3D11_INPUT_ELEMENT_DESC const layoutDesc[],
    _In_ uint32_t layoutDescNumElements,
    _Out_ ID3D11VertexShader** shader,
    _Out_opt_ ID3D11InputLayout** layout
    )
{
    std::vector<byte> bytecode{ m_basicReaderWriter.ReadData(filename) };

    winrt::check_hresult(
        m_d3dDevice->CreateVertexShader(
            bytecode.data(),
            bytecode.size(),
            nullptr,
            shader
            )
        );

    winrt::check_hresult(
        m_d3dDevice->CreateInputLayout(
        layoutDesc,
        layoutDescNumElements,
        bytecode.data(),
        bytecode.size(),
        layout
    ));

    /*
    if (layout != nullptr)
    {
        CreateInputLayout(
            bytecode.data(),
            bytecode.size(),
            layoutDesc,
            layoutDescNumElements,
            layout
            );  
    }
    */
}

winrt::Windows::Foundation::IAsyncAction BasicLoader::LoadShaderAsync(
    _In_ winrt::hstring filename,
    _In_reads_opt_(layoutDescNumElements) D3D11_INPUT_ELEMENT_DESC const layoutDesc[],
    _In_ uint32_t layoutDescNumElements,
    _Out_ ID3D11VertexShader** shader,
    _Out_opt_ ID3D11InputLayout** layout
    )
{
    // This method assumes that the lifetime of input arguments may be shorter
    // than the duration of this operation. In order to ensure accurate results, a
    // copy of all arguments passed by pointer must be made. The method then
    // ensures that the lifetime of the copied data exceeds that of the coroutine.

    // Create copies of the layoutDesc array as well as the SemanticName strings,
    // both of which are pointers to data whose lifetimes may be shorter than that
    // of this method's coroutine.
    std::vector<D3D11_INPUT_ELEMENT_DESC> layoutDescCopy;
    std::vector<std::string> layoutDescSemanticNamesCopy;
    if (layoutDesc != nullptr)
    {
        layoutDescCopy = { layoutDesc, layoutDesc + layoutDescNumElements };
        layoutDescSemanticNamesCopy.reserve(layoutDescNumElements);

        for (auto&& desc : layoutDescCopy)
        {
            desc.SemanticName = layoutDescSemanticNamesCopy.emplace(layoutDescSemanticNamesCopy.end(), desc.SemanticName)->c_str();
        }
    }

    auto bytecode = co_await m_basicReaderWriter.ReadDataAsync(filename);

    winrt::check_hresult(
        m_d3dDevice->CreateVertexShader(
            bytecode.data(),
            bytecode.Length(),
            nullptr,
            shader
            )
        );

    

    if (layout != nullptr)
    {
        CreateInputLayout(
            bytecode.data(),
            bytecode.Length(),
            layoutDesc == nullptr ? nullptr : layoutDescCopy.data(),
            layoutDescNumElements,
            layout
            );

        
    }
}

void BasicLoader::LoadShader(
    _In_ winrt::hstring const& filename,
    _Out_ ID3D11PixelShader** shader
    )
{
    std::vector<byte> bytecode{ m_basicReaderWriter.ReadData(filename) };

    winrt::check_hresult(
        m_d3dDevice->CreatePixelShader(
            bytecode.data(),
            bytecode.size(),
            nullptr,
            shader
            )
        );

    
}

winrt::Windows::Foundation::IAsyncAction BasicLoader::LoadShaderAsync(
    _In_ winrt::hstring filename,
    _Out_ ID3D11PixelShader** shader
    )
{
    auto bytecode = co_await m_basicReaderWriter.ReadDataAsync(filename);

    winrt::check_hresult(
        m_d3dDevice->CreatePixelShader(
            bytecode.data(),
            bytecode.Length(),
            nullptr,
            shader
            )
        );

    
}

void BasicLoader::LoadShader(
    _In_ winrt::hstring const& filename,
    _Out_ ID3D11ComputeShader** shader
    )
{
    std::vector<byte> bytecode{ m_basicReaderWriter.ReadData(filename) };

    winrt::check_hresult(
        m_d3dDevice->CreateComputeShader(
            bytecode.data(),
            bytecode.size(),
            nullptr,
            shader
            )
        );

    
}

winrt::Windows::Foundation::IAsyncAction BasicLoader::LoadShaderAsync(
    _In_ winrt::hstring filename,
    _Out_ ID3D11ComputeShader** shader
    )
{
    auto bytecode = co_await m_basicReaderWriter.ReadDataAsync(filename);

    winrt::check_hresult(
        m_d3dDevice->CreateComputeShader(
            bytecode.data(),
            bytecode.Length(),
            nullptr,
            shader
            )
        );

    
}

void BasicLoader::LoadShader(
    _In_ winrt::hstring const& filename,
    _Out_ ID3D11GeometryShader** shader
    )
{
    std::vector<byte> bytecode{ m_basicReaderWriter.ReadData(filename) };

    winrt::check_hresult(
        m_d3dDevice->CreateGeometryShader(
            bytecode.data(),
            bytecode.size(),
            nullptr,
            shader
            )
        );

    
}

winrt::Windows::Foundation::IAsyncAction BasicLoader::LoadShaderAsync(
    _In_ winrt::hstring filename,
    _Out_ ID3D11GeometryShader** shader
    )
{
    auto bytecode = co_await m_basicReaderWriter.ReadDataAsync(filename);

    winrt::check_hresult(
        m_d3dDevice->CreateGeometryShader(
            bytecode.data(),
            bytecode.Length(),
            nullptr,
            shader
            )
        );

    
}

void BasicLoader::LoadShader(
    _In_ winrt::hstring const& filename,
    _In_reads_opt_(numEntries) const D3D11_SO_DECLARATION_ENTRY* streamOutDeclaration,
    _In_ uint32_t numEntries,
    _In_reads_opt_(numStrides) const uint32_t* bufferStrides,
    _In_ uint32_t numStrides,
    _In_ uint32_t rasterizedStream,
    _Out_ ID3D11GeometryShader** shader
    )
{
    std::vector<byte> bytecode{ m_basicReaderWriter.ReadData(filename) };

    winrt::check_hresult(
        m_d3dDevice->CreateGeometryShaderWithStreamOutput(
            bytecode.data(),
            bytecode.size(),
            streamOutDeclaration,
            numEntries,
            bufferStrides,
            numStrides,
            rasterizedStream,
            nullptr,
            shader
            )
        );

    
}

winrt::Windows::Foundation::IAsyncAction BasicLoader::LoadShaderAsync(
    _In_ winrt::hstring filename,
    _In_reads_opt_(numEntries) const D3D11_SO_DECLARATION_ENTRY* streamOutDeclaration,
    _In_ uint32_t numEntries,
    _In_reads_opt_(numStrides) const uint32_t* bufferStrides,
    _In_ uint32_t numStrides,
    _In_ uint32_t rasterizedStream,
    _Out_ ID3D11GeometryShader** shader
    )
{
    // This method assumes that the lifetime of input arguments may be shorter
    // than the duration of this coroutine. In order to ensure accurate results, a
    // copy of all arguments passed by pointer must be made. The method then
    // ensures that the lifetime of the copied data exceeds that of the coroutine.

    // Create copies of the streamOutDeclaration array as well as the SemanticName
    // strings, both of which are pointers to data whose lifetimes may be shorter
    // than that of this method's coroutine.
    std::vector<D3D11_SO_DECLARATION_ENTRY> streamOutDeclarationCopy;
    std::vector<std::string> streamOutDeclarationSemanticNamesCopy;
    if (streamOutDeclaration != nullptr)
    {
        streamOutDeclarationCopy = { streamOutDeclaration, streamOutDeclaration + numEntries };
        streamOutDeclarationSemanticNamesCopy.reserve(numEntries);

        for (auto&& desc : streamOutDeclarationCopy)
        {
            desc.SemanticName = streamOutDeclarationSemanticNamesCopy.emplace(streamOutDeclarationSemanticNamesCopy.end(), desc.SemanticName)->c_str();
        }
    }

    // Create a copy of the bufferStrides array, which is a pointer to data
    // whose lifetime may be shorter than that of this method's coroutine.
    std::vector<uint32_t> bufferStridesCopy;
    if (bufferStrides != nullptr)
    {
        bufferStridesCopy = { bufferStrides, bufferStrides + numStrides };
    }

    auto bytecode = co_await m_basicReaderWriter.ReadDataAsync(filename);

    winrt::check_hresult(
        m_d3dDevice->CreateGeometryShaderWithStreamOutput(
            bytecode.data(),
            bytecode.Length(),
            streamOutDeclaration == nullptr ? nullptr : streamOutDeclarationCopy.data(),
            numEntries,
            bufferStrides == nullptr ? nullptr : bufferStridesCopy.data(),
            numStrides,
            rasterizedStream,
            nullptr,
            shader
            )
        );

    
}

void BasicLoader::LoadShader(
    _In_ winrt::hstring const& filename,
    _Out_ ID3D11HullShader** shader
    )
{
    std::vector<byte> bytecode{ m_basicReaderWriter.ReadData(filename) };

    winrt::check_hresult(
        m_d3dDevice->CreateHullShader(
            bytecode.data(),
            bytecode.size(),
            nullptr,
            shader
            )
        );

    
}

winrt::Windows::Foundation::IAsyncAction BasicLoader::LoadShaderAsync(
    _In_ winrt::hstring filename,
    _Out_ ID3D11HullShader** shader
    )
{
    auto bytecode = co_await m_basicReaderWriter.ReadDataAsync(filename);

    winrt::check_hresult(
        m_d3dDevice->CreateHullShader(
            bytecode.data(),
            bytecode.Length(),
            nullptr,
            shader
            )
        );

    
}

void BasicLoader::LoadShader(
    _In_ winrt::hstring const& filename,
    _Out_ ID3D11DomainShader** shader
    )
{
    std::vector<byte> bytecode{ m_basicReaderWriter.ReadData(filename) };

    winrt::check_hresult(
        m_d3dDevice->CreateDomainShader(
            bytecode.data(),
            bytecode.size(),
            nullptr,
            shader
            )
        );

    
}

winrt::Windows::Foundation::IAsyncAction BasicLoader::LoadShaderAsync(
    _In_ winrt::hstring filename,
    _Out_ ID3D11DomainShader** shader
    )
{
    auto bytecode = co_await m_basicReaderWriter.ReadDataAsync(filename);

    winrt::check_hresult(
        m_d3dDevice->CreateDomainShader(
            bytecode.data(),
            bytecode.Length(),
            nullptr,
            shader
            )
        );

    
}

