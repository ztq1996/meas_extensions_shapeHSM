#ifndef LSST_MEAS_EXTENSIONS_SHAPEHSM_HSMADAPTER_H
#define LSST_MEAS_EXTENSIONS_SHAPEHSM_HSMADAPTER_H

#include "lsst/afw/image/Image.h"
#include "galsim/Image.h"

namespace lsst { namespace meas { namespace extensions { namespace shapeHSM {

/// Convert from afw Image to GalSim's ImageView
template <typename PixelT>
class ImageConverter
{
public:
    /// Ctor
    ///
    /// ImageView wants us to provide an "owner" for the pixels (so it can hold a
    /// reference to keep them from being destroyed underneath it?), but ndarray
    /// doesn't allow us to access a shared_ptr to the pixels.  Instead, we use a
    /// dummy shared ptr to a single (new) pixel.  The ImageConverter is holding
    /// on to this too (RAII) so we shouldn't have any memory worries.
    ImageConverter(PTR(afw::image::Image<PixelT>) image, afw::geom::Box2I box) :
        _image(image), _owner(new PixelT), _box(box) {}
    ImageConverter(PTR(afw::image::Image<PixelT>) image) :
        _image(image), _owner(new PixelT), _box(image->getBBox(afw::image::LOCAL)) {}

    /// Conversion
    galsim::ImageView<PixelT> getImageView() const {
        // Make upper bound inclusive: that's the GalSim standard, but not LSST
        galsim::Bounds<int> const bounds(_box.getMinX(), _box.getMaxX() - 1,
                                         _box.getMinY(), _box.getMaxY() - 1);
        typename afw::image::Image<PixelT>::Array array = _image->getArray();
        int const stride = array.template getStride<0>();
        // Advance into the array by the box, as GalSim's BaseImage.addressPixel takes it off
        // Remove xy0, as that's just on there for convenience (e.g., tract/patch), or it's included
        // in the box and stride.
        PixelT* ptr = reinterpret_cast<PixelT*>(array.getData() + (_box.getMinY() - _image->getY0())*stride +
                                                _box.getMinX() - _image->getX0());
        return galsim::ImageView<PixelT>(ptr, _owner, stride, bounds);
    }

#if 0
    /// Write to FITS
    ///
    /// Intended for debugging, as it inspects the ImageView
    void writeFits(std::string const& filename) const {
        galsim::ImageView<PixelT> view = getImageView();
        galsim::Bounds<int> const& bounds = view.getBounds();
        int const width = view.getXMax() - view.getXMin() + 1;
        int const height = view.getYMax() - view.getYMin() + 1;
        afw::image::Image<PixelT> out(width, height);
        int const x0 = bounds.getXMin(), y0 = bounds.getYMin();
        for (int y = bounds.getYMin(); y <= bounds.getYMax(); ++y) {
            for (int x = bounds.getXMin(); x <= bounds.getXMax(); ++x) {
                out(x - x0, y - y0) = view.at(x, y);
            }
        }
        out.writeFits(filename);
    }
    void writeFits(boost::format const& format) const { writeFits(format.str()); }
#endif

private:
    PTR(afw::image::Image<PixelT>) _image;
    boost::shared_ptr<PixelT> _owner;
    afw::geom::Box2I _box;
};


/// Convert an afw Mask into an image for use with HSM
///
/// HSM uses a mask where 0 = bad, 1 = good.
inline
PTR(afw::image::Image<int>) convertMask(
    afw::image::Mask<afw::image::MaskPixel> const& afwMask, ///< Traditional afw Mask using mask planes
    afw::geom::Box2I const& bbox,            ///< Bounding box of interest
    afw::image::MaskPixel const badPixelMask ///< Mask for selecting bad pixels
    )
{
    typedef afw::image::Image<int> ImageI;
    typedef afw::image::Mask<afw::image::MaskPixel> Mask;
    PTR(ImageI) hsmMask = boost::make_shared<ImageI>(bbox.getDimensions());
    for (int y = 0; y < bbox.getHeight(); ++y) {
        Mask::const_x_iterator in = afwMask.x_at(bbox.getMinX() - afwMask.getX0(), y);
        ImageI::x_iterator out = hsmMask->row_begin(y);
        ImageI::const_x_iterator end = hsmMask->row_end(y);
        for (; out != end; ++in, ++out) {
            *out = (*in & badPixelMask) ? 0 : 1;
        }
    }
    return hsmMask;
}

}}}} // namespace lsst::meas::extensions::shapeHSM


#endif // ifndef LSST_MEAS_EXTENSIONS_SHAPEHSM_HSMADAPTER_H
