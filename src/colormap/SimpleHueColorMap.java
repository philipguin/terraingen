package colormap;
import java.awt.Color;



public class SimpleHueColorMap implements Mapper<Float, Color>
{
	private final float startHue, endHue, saturation, brightness;
	
	public SimpleHueColorMap(float startHue, float endHue, float saturation, float brightness)
	{
		this.startHue = startHue;
		this.endHue = endHue;
		this.saturation = saturation;
		this.brightness = brightness;
	}
	
	@Override
	public Color map(Float value)
	{
		return new Color(Color.HSBtoRGB(startHue + (endHue - startHue) * value, saturation, brightness));
	}
}
